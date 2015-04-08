# -*- coding: utf-8 -*-

## creating a temporary reStructeredTxt version of README for setup.py to read in
def convert_md2rst():
   try:
      import pypandoc
      description_rst = pypandoc.convert('README.md', 'rst')
      with open('README.rst', 'a') as f:
         f.write(description_rst)
   except (IOError, ImportError):
      description_rst = open('README.md').read()

def main():
   convert_md2rst()

if __name__ == '__main__':
    main()
