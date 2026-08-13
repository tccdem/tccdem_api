"""
tccdem_api : API for TCCDEM Database
"""

# 1. Package Metadata
__version__ = "0.1.0"
__author__ = "Cem Oran"

# 2. Expose core functions/classes at the root package level
from .tccdem_api import MaterialDatabaseAPI
# from PNcsp_Plus.db import DBconnector 


# 3. Control public exports for wildcard imports (`from PNcsp_Plus import *`)
__all__ = ["MaterialDatabaseAPI"]