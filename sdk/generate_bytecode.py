import os
import compileall

assert os.path.exists('setup.py') and os.path.exists('ComPASS') and os.path.isdir('ComPASS')

compileall.compile_dir('ComPASS')

