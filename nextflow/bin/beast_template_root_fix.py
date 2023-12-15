#!/usr/bin/env python
# import click
import argparse

parser = argparse.ArgumentParser(description='A script that takes in string and file arguments')
parser.add_argument('--templated_beast_xml_path', type=str, help='')
parser.add_argument('--out_path', type=str, help='')




#@click.command()
#@click.argument("templated_beast_xml_path")
#@click.argument("out_path")
def run(templated_beast_xml_path, out_path):
    """
    Read in a templated beast xml file that is missing the '<taxa id="ingroup">'
    tag (used to force rooting on the naive seqeunce), and write out the xml
    file with the tag added.

    Parameters:
        templated_beast_xml_path (str): The path of the templated beast xml file. The
            file should be fully templated and formatted to root on the naive sequence,
            except that it is missing the xml tag defining the ingroup.
        out_path (str): The path to write the fully formatted beast xml file.
    """
    with open(templated_beast_xml_path) as the_file:
        lines = the_file.readlines()

    start_tag = '<taxa id="taxa">'
    end_tag = "</taxa>"
    start = next(j for j, line in enumerate(lines) if line.find(start_tag) != -1)
    stop = next(j for j, line in enumerate(lines) if line.find(end_tag) != -1)

    new_start_tag = lines[start].replace('"taxa"', '"ingroup"')
    new_end_tag = lines[stop]

    to_keep, to_skip = "taxon id", 'id="naive@0"'
    use_line = lambda s: s.find(to_keep) != -1 and s.find(to_skip) == -1
    edit_line = lambda s: s.replace("id", "idref", 1).replace(">", "/>", 1)
    lines_to_add = (
        edit_line(line) for line in lines[start + 1 : stop] if use_line(line)
    )

    with open(out_path, "w") as the_file:
        for line in lines[: stop + 1]:
            the_file.write(line)
        the_file.write(new_start_tag)
        for line in lines_to_add:
            the_file.write(line)
        the_file.write(new_end_tag)
        for line in lines[stop + 1 :]:
            the_file.write(line)

    return None


if __name__ == "__main__":
    
    args = parser.parse_args()

    run(args.templated_beast_xml_path, args.out_path)
    





