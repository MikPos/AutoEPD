import mod
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import ap
from test_book import ExM

ExM.run(prettyPrint=True)