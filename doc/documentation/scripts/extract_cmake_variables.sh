#!/bin/bash

# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# This script extracts CMake variables from the CMakeCache.txt file.
CACHE_FILE=$1

awk '
BEGIN {
  print "Variable,Type,Description"
  doc = ""
}
/^\/\// {
  line = $0
  gsub("^//", "", line)
  gsub(/^[ \t]+|[ \t]+$/, "", line)
  doc = (doc == "") ? line : doc " " line
  next
}
/^[^#].*:/ {
  split($0, a, ":")
  name = a[1]
  split(a[2], t, "=")
  type = t[1]
  if (name ~ /^FOUR_C/ && type != "INTERNAL") {
    gsub(/"/, "\"\"", doc)
    printf "\"``%s``\",\"%s\",\"%s\"\n", name, type, doc
  }
  doc = ""  # reset for next block
}
' ${CACHE_FILE}
