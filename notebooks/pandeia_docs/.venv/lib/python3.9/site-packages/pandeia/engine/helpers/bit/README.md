# TEMPORARY MIGRATION CODE !!!

The vast majority of this directory is TEMPORARY for the conversion of APT queries.

This is NOT run during normal engine calculations and once we have developed a pandeia-native tool, most of this content will go away.

When this folder is cleaned out, keep any Python scripts with names matching "bit*.py"

History:
============================================================

* 20211020: This directory was created simply to facilitate running APT BIT queries through the pandeia engine.  A smaller share of this dir is needed inoder to go from APT query->pweb.  A larger chunk is needed for the pweb->peng conversion.


* Needed for just query->pweb:

```
__init__.py
bit.py
compute.py
config/
input_conversion.py
pyetc_form_defaults/
pyetc_util.py
translate.py
```

* Needed for pweb->peng:

```
etc_web_data/
extraction.py
instruments/ # this is a LOT of content!
perform_calculation.py
sky.py
source.py
string_constants.py
and maybe more ...
```
