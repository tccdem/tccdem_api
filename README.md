# Material Database API

A Python interface for storing, querying, and managing crystal structure data and calculated material properties.

---

## Features

- **Compressed Storage:** Automatically serializes and compresses `pymatgen` `Structure` objects (`zlib`) to keep storage lightweight.
- **Uniqueness Checking:** Compares candidate structures against existing records in the database using `StructureMatcher`.
- **Property Tracking:** Associates Machine Learning properties  with stored structures.
- **CIF Export:** Converts stored structures back into symmetrized CIF files.
- **Flexibility:** Query by composition and project name.

---

## Configuration & Environment Setup

The API automatically loads environment variables from your home directory at `~/.config/tccdem_db/db.env`. Ensure your `.env` file contains the following keys:

```env
tccdem_db_host=localhost
tccdem_db_user=your_username
tccdem_db_pswd=your_password
tccdem_db_name=your_db_name
tccdem_db_port=3306
```

## Usage
### Initialization
```python
from tccdem_db import MaterialDatabaseAPI

# Initialize using environment variables
db = MaterialDatabaseAPI()
```
### Get Data
Get structure data :
```python
# Get all data
structures = db.get_structure()

# Get data by filtering (by -> {'composition','project','user'})
structures = db.get_structure(by=str, entry=str)
```

Get property data by filtering:
```python
# Get all data
properties = db.get_property()

# Get data by filtering (by -> {'composition','project','user'})
properties = db.get_property(by=str, entry=str)
```
Get user and project names:
```python
projects = db.get_projects()
users = db.get_users()
```
Fetch all rows in tables:
```python
# table_name -> {'Uploads', 'Compositions', 'Structures', 'Properties'}
table = db.fetch_table(table_name=str)
```
Write CIF to a specific path:
```python
file_path = db.write_cif(struc_id=int, path=str)
```
Generate Castep input in a specific path:
```python
db.generate_castep_input(struc_id=int, path=str)
```
### Upload Data
Upload structure with property (if no property, leave properties_dict as an empty dict) :
```python
struc_id = db.upload_data(
    project=str,
    structure_obj=pmg_structure_object,
    properties_dict={
    "property_name1": {"value":float, "unit":str, "program":str},
    "property_name2": {"value":float, "unit":str, "program":str}
    },
    dimension=int,
    prototype=str,
    is_struct_new=bool,
    is_sg_new=bool,
    relaxation=str
)
```

Upload property for an existing structure:
```python
new_upload_id = db.upload_property(
    struc_id=int,
    properties_dict={
    "property_name1": {"value":float, "unit":str, "program":str},
    "property_name2": {"value":float, "unit":str, "program":str}
    },
    project=str
)
```

### Remove Data
```python
# Remove by filterin(by -> {'struc_id', 'prop_id', 'project','user'})
db.clear_entries(by=str, entry=int)

# Remove all data in the database
db.clear_database()
```
