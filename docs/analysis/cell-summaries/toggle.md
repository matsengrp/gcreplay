# Cell Summaries with Toggle

<div class="toggle-buttons">
  <button class="notebook-toggle-btn active" data-target="default-notebook">Default</button>
  <button class="notebook-toggle-btn" data-target="naive-first-notebook">Naive Reversions First</button>
  <button class="notebook-toggle-btn" data-target="naive-no-bp-notebook">Naive Reversions No BP</button>
</div>

<style>
.toggle-buttons {
  margin-bottom: 20px;
}
.notebook-toggle-btn {
  padding: 8px 16px;
  margin-right: 8px;
  border: 1px solid #ddd;
  border-radius: 4px;
  background: #f8f8f8;
  cursor: pointer;
}
.notebook-toggle-btn.active {
  background: #2196F3;
  color: white;
  border-color: #2196F3;
}
.notebook-variation {
  display: none;
}
</style>

<div id="default-notebook" class="notebook-variation">
  <iframe src="default.ipynb" width="100%" height="800px" frameborder="0"></iframe>
</div>

<div id="naive-first-notebook" class="notebook-variation">
  <iframe src="../analysis/cell-summaries/naive_reversions_first.ipynb" width="100%" height="800px" frameborder="0"></iframe>
</div>

<div id="naive-no-bp-notebook" class="notebook-variation">
  <iframe src="../analysis/cell-summaries/naive_reversions_no_bp.ipynb" width="100%" height="800px" frameborder="0"></iframe>
</div>
