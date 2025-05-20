import stat
import numpy as np
import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
import equinox as eqx
import abc
import diffrax as dx

from typing import Callable, Self
from jaxtyping import Array, Float, Scalar, PyTree

import matplotlib.pyplot as plt


# NB: It's not possible to have static JAX arrays. That's a requirement from JAX.


class AbstractSpatialDiscretization(eqx.Module):
    x0: float = eqx.field(static=True)
    x_final: float = eqx.field(static=True)
    n: int = eqx.field(static=True)
    vals: eqx.AbstractVar[Float[Array, "n"]]

    @property
    def δx(self):
        return (self.x_final - self.x0) / (self.n - 1)

    @property
    def x(self):
        return jnp.linspace(self.x0, self.x_final, self.n)

    def check_aligned(self, other: Self) -> None:
        r"""Check if another field is spatially aligned with this one (i.e., it has the
        same domain parameters).

        Args:
            other: Other field.

        Raises:
            eqx.EquinoxTracetimeError: If ``other`` is not aligned.
        """
        eqx.error_if(
            other.vals,
            jnp.logical_not(jnp.array_equal(self.x0, other.x0)),
            f"Mismatched lower bounds {self.x0} and {other.x0}",
        )
        eqx.error_if(
            other.vals,
            jnp.logical_not(jnp.array_equal(self.x_final, other.x_final)),
            f"Mismatched upper bounds {self.x_final} and {other.x_final}",
        )
        eqx.error_if(
            other.vals,
            jnp.logical_not(jnp.array_equal(self.n, other.n)),
            f"Mismatched number of points {self.n} and {other.n}",
        )

    def map(self, fn):
        return SpatialDiscretization(self.x0, self.x_final, self.n, fn(self.vals))

    def binop(self, other, fn):
        if isinstance(other, AbstractSpatialDiscretization):
            self.check_aligned(other)
            other = other.vals
        return SpatialDiscretization(
            self.x0, self.x_final, self.n, fn(self.vals, other)
        )

    def integral(self):
        return jax.scipy.integrate.trapezoid(self.vals, dx=self.δx, axis=-1)

    def indefinite_integral(self):
        return self.map(
            lambda vals: jnp.concatenate(
                [
                    jnp.zeros_like(vals[..., :1]),
                    jnp.cumsum(self.δx * (vals[..., :-1] + vals[..., 1:]) / 2, axis=-1),
                ],
                axis=-1,
            )
        )

    def convolve(self, other):
        if isinstance(other, AbstractSpatialDiscretization):
            self.check_aligned(other)
            other = other.vals
        return SpatialDiscretization(
            self.x0,
            self.x_final,
            self.n,
            jnp.convolve(self.vals, other, mode="same") * self.δx,
        )

    def __add__(self, other):
        return self.binop(other, lambda x, y: x + y)

    def __sub__(self, other):
        return self.binop(other, lambda x, y: x - y)

    def __mul__(self, other):
        return self.binop(other, lambda x, y: x * y)

    def __truediv__(self, other):
        return self.binop(other, lambda x, y: x / y)

    def __radd__(self, other):
        return self.binop(other, lambda x, y: y + x)

    def __rsub__(self, other):
        return self.binop(other, lambda x, y: y - x)

    def __rmul__(self, other):
        return self.binop(other, lambda x, y: y * x)

    def __rtruediv__(self, other):
        return self.binop(other, lambda x, y: y / x)

    def plot(self, axis=None, t=None, **kwargs):
        if axis is not None:
            plt.sca(axis)
        if self.vals.ndim == 2:
            return plt.pcolor(t, self.x, self.vals.T, **kwargs)
        return plt.plot(self.x, self.vals, **kwargs)


class SpatialDiscretization(AbstractSpatialDiscretization):
    vals: Float[Array, "n"] = eqx.field(converter=jnp.asarray)

    def __check_init__(self):
        eqx.error_if(
            self.vals,
            jnp.logical_not(jnp.array_equal(self.vals.shape[-1], self.n)),
            f"Domain has size {self.vals.shape[-1]} but should have size {self.n}",
        )

    def eval(self, x: Float[Array, "n"]) -> Float[Array, "n"]:
        return jnp.interp(x, self.x, self.vals)


class ParameterizedSpatialDiscretization(AbstractSpatialDiscretization, abc.ABC):
    @property
    def vals(self):
        return self.eval(self.x)

    @abc.abstractmethod
    def eval(self, x: Float[Array, "n"]) -> Float[Array, "n"]:
        pass

    def __repr__(self):
        return eqx.tree_pformat(
            self,
            short_arrays=False,
            truncate_leaf=lambda x: isinstance(x, jnp.ndarray) and x.size > 100,
        )


class Sigmoid(ParameterizedSpatialDiscretization):
    α: float = eqx.field(converter=jnp.asarray)
    β: float = eqx.field(converter=jnp.asarray)
    γ: float = eqx.field(converter=jnp.asarray)
    ν: float = eqx.field(converter=jnp.asarray)

    def eval(self, x):
        # return self.α * (jnp.tanh(self.β * x + self.γ) - jnp.tanh(self.γ))
        return self.α / (1 + jnp.exp(-(self.β * x - self.γ)))**self.ν

    def derivative(self):
        return SpatialDiscretization(
            self.x0,
            self.x_final,
            self.n,
            self.α * self.β * self.ν * jnp.exp(self.β * self.ν * self.x + self.γ) * (jnp.exp(self.γ) + jnp.exp(self.β * self.x)) ** (-self.ν - 1)
        )


class Gaussian(ParameterizedSpatialDiscretization):
    mass: float = eqx.field(converter=jnp.asarray)
    mean: float = eqx.field(converter=jnp.asarray)
    std: float = eqx.field(converter=jnp.asarray)

    def __check_init__(self):
        eqx.error_if(
            self.std,
            jnp.any(self.std <= 0),
            f"Standard deviations must be positive",
        )

    def eval(self, x):
        return (
            self.mass
            * jnp.exp(-0.5 * ((x - self.mean) / self.std) ** 2)
            / (self.std * jnp.sqrt(2 * jnp.pi))
        )
    
    def logpdf(self, x):
        return -0.5 * ((x - self.mean) / self.std) ** 2 - jnp.log(self.std * jnp.sqrt(2 * jnp.pi))

    def integral(self):
        return self.mass


class GaussianMixture(ParameterizedSpatialDiscretization):
    weights: Float[Array, "r"] = eqx.field(converter=jnp.asarray)
    means: Float[Array, "r"] = eqx.field(converter=jnp.asarray)
    stds: Float[Array, "r"] = eqx.field(converter=jnp.asarray)

    def __check_init__(self):
        eqx.error_if(
            self.stds,
            jnp.any(self.stds <= 0),
            f"Standard deviations must be positive",
        )

    def eval(self, x):
        return jnp.sum(
            self.weights
            * jnp.exp(-0.5 * ((x[:, None] - self.means) / self.stds) ** 2)
            / (self.stds * jnp.sqrt(2 * jnp.pi)),
            axis=-1,
        )
    
    def logpdf(self, x):
        weights = self.weights / jnp.sum(self.weights)
        return jax.scipy.special.logsumexp(
            jnp.log(weights)
            - 0.5 * ((x[:, None] - self.means) / self.stds) ** 2
            - jnp.log(self.stds * jnp.sqrt(2 * jnp.pi)),
            axis=-1,
        )

    def integral(self):
        return jnp.sum(self.weights)

class Interp(ParameterizedSpatialDiscretization):
    xp: Float[Array, "k"] = eqx.field(static=True, converter=jnp.asarray)
    yp: Float[Array, "k"] = eqx.field(converter=jnp.asarray)

    def eval(self, x):
        return jnp.interp(x, self.xp, self.yp)


def _is_scalar_field(node):
    return isinstance(node, SpatialDiscretization)


def _field_structure(field):
    return jax.tree_util.tree_structure(field, is_leaf=_is_scalar_field)


class PDETerm(dx.AbstractTerm):
    r"""A term representing :math:`f(t, y(x, t), args) \mathrm{d}t`. That is to say,
    the term appearing on the right hand side of a PDE, in which the control is time.

    ``vector_field`` should return some PyTree, with the same structure as the initial
    state ``y0``, and with every leaf broadcastable to the equivalent leaf in ``y0``.
    """

    vector_field: Callable[
        [Scalar, SpatialDiscretization, PyTree], SpatialDiscretization
    ]

    def vf(
        self, t: Scalar, y: SpatialDiscretization, args: PyTree
    ) -> SpatialDiscretization:
        out = self.vector_field(t, y, args)
        if _field_structure(out) != _field_structure(y):
            raise ValueError(
                f"Vector field output structure {_field_structure(out)} "
                f"does not match input structure {_field_structure(y)}"
            )
        return out
        # return jax.tree_util.tree_map(lambda o, yi: jnp.broadcast_to(o, jnp.shape(yi)), out, y)

    @staticmethod
    def contr(t0: Scalar, t1: Scalar) -> Scalar:
        return t1 - t0

    @staticmethod
    def prod(vf: SpatialDiscretization, control: Scalar) -> PyTree:
        return jax.tree_util.tree_map(lambda v: control * v, vf)


class CrankNicolson(dx.AbstractSolver):
    rtol: float
    atol: float

    term_structure = PDETerm
    interpolation_cls = dx.LocalLinearInterpolation

    def order(self, terms):
        return 2

    def init(self, terms, t0, t1, y0, args):
        return None

    def step(self, terms, t0, t1, y0, args, solver_state, made_jump):
        del solver_state, made_jump
        δt = t1 - t0
        f0 = terms.vf(t0, y0, args)

        def keep_iterating(val):
            _, not_converged = val
            return not_converged

        def fixed_point_iteration(val):
            y1, _ = val
            new_y1 = y0 + 0.5 * δt * (f0 + terms.vf(t1, y1, args))
            diff = jnp.abs((new_y1 - y1).vals)
            max_y1 = jnp.maximum(jnp.abs(y1.vals), jnp.abs(new_y1.vals))
            scale = self.atol + self.rtol * max_y1
            not_converged = jnp.any(diff > scale)
            return new_y1, not_converged

        euler_y1 = y0 + δt * f0
        y1, _ = eqx.internal.while_loop(
            keep_iterating,
            fixed_point_iteration,
            (euler_y1, False),
            kind="checkpointed",
            checkpoints=100,
        )

        y_error = y1 - euler_y1
        dense_info = dict(y0=y0, y1=y1)

        solver_state = None
        result = dx.RESULTS.successful
        return y1, y_error, dense_info, solver_state, result

    def func(self, terms, t0, y0, args):
        return terms.vf(t0, y0, args)
