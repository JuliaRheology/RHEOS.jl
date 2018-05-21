#

# for several months, at least, emacs regularly insersts whitespace at the end of lines

find . -type f -name "*.jl README* REQU* .travis.yml LICENSE* " -exec sh -c 'for i;do sed 's/[[:space:]]*$//' "$i">/tmp/.$$ && cat /tmp/.$$ > "$i";done' arg0 {} +
