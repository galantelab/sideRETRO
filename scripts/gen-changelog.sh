#!/bin/sh
#
# sideRETRO - A pipeline for detecting Somatic Insertion of DE novo RETROcopies
# Copyright (C) 2019-2026 Thiago L. A. Miller <tmiller@mochsl.org.br
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Generate a GNU style ChangeLog from the git history.
#
# Usage: gen-changelog.sh [SOURCE_ROOT] [OUTPUT]
#
# Without arguments, run as a 'meson dist' script: Read the
# repository at $MESON_SOURCE_ROOT and write $MESON_DIST_ROOT/ChangeLog

set -e

source_root="${1:-$MESON_SOURCE_ROOT}"
output="${2:-$MESON_DIST_ROOT/ChangeLog}"

if [ -z "$source_root" ] || [ -z "$output" ]; then
	echo "Usage: $0 [SOURCE_ROOT] [OUTPUT]" >&2
	exit 1
fi

git -C "$source_root" log \
	--date=short \
	--pretty=format:"%ad  %an  <%ae>%n%n%w(76,8,8)%B" \
	> "$output"
