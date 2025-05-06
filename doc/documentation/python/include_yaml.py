# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

import yaml
import pathlib
import shutil
import jinja2
import argparse

PATH_TO_TESTS = ""
RECOGNIZED_TYPES = ["rst", "md"]


def load_input_file(yaml_file):
    data = yaml.safe_load((PATH_TO_TESTS / pathlib.Path(yaml_file)).read_text())
    return data


def yaml_dump(data, filetype="rst"):
    if filetype == "rst":
        rststring = ""
        yaml_data = yaml.safe_dump(data, sort_keys=False).split("\n")
        for line in yaml_data:
            rststring += "    " + line + "\n"
        rststring += "\n"
        return ".. code-block:: yaml\n\n" + rststring + "\n"
    elif filetype == "md":
        return "```yaml\n" + yaml.safe_dump(data, sort_keys=False) + "```\n"
    else:
        raise TypeError(f"Filetype {filetype} for yaml_dump cannot be recognized yet.")


def section_dump(input_file_data, section_names, filetype="rst"):
    # section_dump can take either one section name or a list of section names
    if isinstance(section_names, str):
        section_names = [section_names]
    yaml_dict = {}
    for section in section_names:
        yaml_dict[section] = input_file_data[section]
    return yaml_dump(yaml_dict, filetype)


def convert(template_path, rendering_path, input_file_path):
    global PATH_TO_TESTS
    PATH_TO_TESTS = input_file_path
    target_dir = pathlib.Path(rendering_path)

    for template_file in pathlib.Path(template_path).glob("*.j2"):

        template = jinja2.Template(template_file.read_text())
        tutorial_name = template_file.stem
        tutorial_rst = target_dir / tutorial_name
        print(f"source: {tutorial_name}, target: {tutorial_rst}")
        tutorial_rst.write_text(
            template.render(
                section_dump=section_dump,
                load_input_file=load_input_file,
                yaml_dump=yaml_dump,
                len=len,
            )
        )
    for filetype in RECOGNIZED_TYPES:
        for file in pathlib.Path(template_path).glob("*." + filetype):
            shutil.copy(file, target_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render all tutorials in the given directory."
    )
    parser.add_argument(
        "template_path", type=str, help="Path to the tutorial templates"
    )
    parser.add_argument(
        "rendering_path", type=str, help="Path to the final tutorial files"
    )
    parser.add_argument("input_file_path", type=str, help="Path to the input files")
    # Parse the arguments
    args = parser.parse_args()

    convert(args.template_path, args.rendering_path, args.input_file_path)
