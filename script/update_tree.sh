
#!/usr/bin/env bash
set -euo pipefail
README="README.md"
TREE=$(find . -maxdepth 3 \
  -path './.git' -prune -o \
  -path './data raw' -prune -o \
  -path './data clean' -prune -o \
  -print | sed 's|^\./||' | sort)
FORMAT=$(echo "$TREE" | awk -F'/' '
{
  indent=""
  for (i=1; i<NF; i++) indent=indent"  "
  print indent"- "$NF
}')
if [ ! -f "$README" ]; then
  cat > "$README" <<'MD'
# ecoChina2

## Repo structure (auto-generated)

<!-- TREE:START -->
<!-- TREE:END -->
MD
fi
perl -0777 -i -pe '
s/<!-- TREE:START -->.*?<!-- TREE:END -->/<!-- TREE:START -->\n```text\n'"$FORMAT"'\n```\n<!-- TREE:END -->/s
' "$README"
EOF

chmod +x script/update_tree.sh
