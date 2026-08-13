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