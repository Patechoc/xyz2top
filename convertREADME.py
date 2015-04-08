# -*- coding: utf-8 -*-

## creating a temporary reStructeredTxt version of README for setup.py to read in
def convert_md2rst(readmeMD="README.md",readmeRST="README.rst"):
   description_rst = ''
   try:
      import pypandoc
      description_rst = pypandoc.convert(readmeMD, 'rst')
      with open(readmeRST, 'w') as f:
         f.write(description_rst)
   except (IOError, ImportError):
      description_rst = open('README.md').read()
   return description_rst

def convert_rst2md(readmeRST="README.rst",readmeMD="README.md"):
   description_md = ''
   try:
      import pypandoc
      description_md = pypandoc.convert(readmeRST, 'md')
      with open(readmeMD, 'w') as f:
         f.write(description_md)
   except (IOError, ImportError):
      description_md = open('README.rst').read()
   return description_md

def main():
   convert_md2rst()
   
if __name__ == '__main__':
    main()
