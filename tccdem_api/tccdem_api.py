import json
import zlib
import os
from pathlib import Path
import mysql.connector
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def load_config():
    from dotenv import load_dotenv
    load_dotenv(Path.home() / ".config/tccdem_db/db.env")

load_config()
PACKAGE_DIR = Path(__file__).resolve().parent


class MaterialDatabaseAPI:
    def __init__(
        self, 
        host: str = os.environ.get("tccdem_db_host"), 
        user: str = os.environ.get("tccdem_db_user"),
        password: str = os.environ.get("tccdem_db_pswd"),
        name: str = os.environ.get("tccdem_db_name"),
        port: str = os.environ.get("tccdem_db_port"),
        use_pure: bool = True,
        ):
        self.config = {
            'host': host,
            'user': user,
            'password': password,
            'database': name,
            'host': host,
            'port': port,
            'charset': 'utf8mb4'
        }

    def _get_connection(self):
        return mysql.connector.connect(**self.config)
    
    def order_formula(self, comp):
        import re

        parts = re.findall(r"([A-Z][a-z]?)(\d*)", comp.strip())
        comp_w1=""
        for el, num in parts:
            comp_w1+=el
            comp_w1+=num if num else "1"

        comp_decomp=[x for x in re.split(r'(?=[A-Z])', comp_w1) if x]
        comp_decomp=sorted(comp_decomp)
        comp_ordered=""
        for j in range(len(comp_decomp)):
            comp_ordered+=comp_decomp[j]
        return comp_ordered        
    
    def StructureControl(self,comp,new_struc,project):
        from pymatgen.analysis.structure_matcher import StructureMatcher
        from pymatgen.core import Structure
        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        query = """
            SELECT 
                s.struc_id,
                s.structure AS compressed_structure
            FROM Structures s
            JOIN Compositions c ON s.comp_id = c.comp_id
            JOIN Uploads u ON u.upload_id = s.upload_id
            WHERE c.formula = %s and u.project = %s;
            """

        try:
            # print("\nStructure Control Module is activated...")
            comp=self.order_formula(comp)
            cursor.execute(query,(comp,project,))
            rows = cursor.fetchall()
            data={}
            for row in rows:
                # 1. Decompress and deserialize the Pymatgen structure
                decompressed_bytes = zlib.decompress(row['compressed_structure'])
                structure_dict = json.loads(decompressed_bytes.decode('utf-8'))
                struct_obj = Structure.from_dict(structure_dict)
                data[row["struc_id"]]=struct_obj

            for key in data.keys():
                sm = StructureMatcher()
                if sm.fit(new_struc, data[key]):
                    similar_to=key
                    print("The structure is similar to struc_id:",key,"\n")
                    return(0)
            # print("The structure is unique.")
            return(1)
                     

        except mysql.connector.Error as err:
                    print(f"Database query failed. Error: {err}")
                    raise err
        finally:
            cursor.close()
            conn.close()
        return data

    def upload_data(self, project, structure_obj: Structure, 
                           properties_dict: dict, dimension: int, prototype: str = None,
                           is_struct_new: bool = None, is_sg_new: bool = None, 
                           similar_to: str = None,relaxation: str = None):
        """
        Uploads an entire CSP generation run to the database in a safe transaction.

        structure_obj expect Ase Structure object

        properties_dict expects format (If there is no property, leave it None): 
        {
            'property_name': {'value': float, 'unit': str, 'program': str}, ...
        }
        """
        def get_PN(Symb):
            csv_path = PACKAGE_DIR / "data" / "Z_PN_Elem_extended_MD.csv"
            with open(csv_path,'r', encoding='utf-8-sig') as file:
                Comp_list={}
                for line in file:
                    if 'PN' in line:
                        continue
                    info=line.strip().replace('\ufeff','').split(',')
                    Comp_list[info[1]]=int(info[0])
            return(Comp_list[Symb])
        
        username=os.environ.get("USER")
        
        conn = self._get_connection()
        cursor = conn.cursor()
        
        try:
            # 1. Prepare Data
            # 1.1 Extract and Handle Composition via Pymatgen
            comp = structure_obj.composition
            formula = self.order_formula(comp.reduced_formula)      ## make it alphabetical

            # Uniqueness Check
            if not self.StructureControl(formula,structure_obj,project):
                return None

            generic = comp.anonymized_formula
            elements = sorted([str(el) for el in comp.elements])
            pn_list = [get_PN(el) for el in elements] ## Open the Z table
            ntypes = len(elements)

            sga = SpacegroupAnalyzer(structure_obj, symprec=1e-2, angle_tolerance=5)
            spacegroup=sga.get_space_group_number()

            natoms=int(structure_obj.num_sites)
            fu=structure_obj.composition.get_reduced_composition_and_factor()[1]

            # 1.2 Serialize and compress Pymatgen structure to Binary MediumBlob
            structure_dict = structure_obj.as_dict()
            json_bytes = json.dumps(structure_dict).encode('utf-8')
            compressed_structure = zlib.compress(json_bytes)

            # 2. Insert into Uploads Table
            upload_query = "INSERT INTO Uploads (username, project) VALUES (%s, %s);"
            cursor.execute(upload_query, (username, project))
            upload_id = cursor.lastrowid

            # 3. Insert into Compounds Table
            comp_query = """
            INSERT INTO Compositions (formula, generic, ntypes, element_list, PN_list)
            VALUES (%s, %s, %s, %s, %s)
            ON DUPLICATE KEY UPDATE comp_id=LAST_INSERT_ID(comp_id);
            """
            cursor.execute(comp_query, (formula, generic, ntypes, json.dumps(elements), json.dumps(pn_list)))
            comp_id = cursor.lastrowid

            # 4. Insert into Structures Table
            struc_query = """
            INSERT INTO Structures (
                comp_id, upload_id, spacegroup, prototype, natoms, fu, 
                dimension, is_struct_new, is_SG_new, similar_to, structure, relaxation
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
            """
            struc_data = (
                comp_id, upload_id, spacegroup, prototype, natoms, 
                fu, dimension, is_struct_new, is_sg_new, similar_to, compressed_structure, relaxation
            )
            cursor.execute(struc_query, struc_data)
            struc_id = cursor.lastrowid

            # 5. Insert into Porperties Table (e.g., MACE Energy)
            if properties_dict:
                prop_query = """
                INSERT INTO Properties (struc_id, upload_id, property, value, unit, program)
                VALUES (%s, %s, %s, %s, %s, %s);
                """
                prop_payload = []
                for prop_name, data in properties_dict.items():
                    prop_payload.append((
                        struc_id, upload_id, prop_name, 
                        data['value'], data.get('unit'), data.get('program')
                    ))
                cursor.executemany(prop_query, prop_payload)

            # Commit if everything succeeds
            conn.commit()
            print(f"Successfully uploaded structure (ID: {struc_id}) and properties.")
            return struc_id

        except mysql.connector.Error as err:
            conn.rollback() # The database is restored to its exact state before the transaction began, preventing partial, corrupted, or inconsistent data.
            print(f"Transaction failed, rolled back changes. Error: {err}")
        except TypeError:
             print("Error: Pymatgen's get_space_group_number() function failed: 'NoneType' object is not subscriptable\n")
             conn.rollback()
        except Exception as e:
            # Catches broader symmetry analysis errors (e.g., spglib / tolerance errors)
            print(f"Error: {e}\n")
            conn.rollback()
        finally:
            cursor.close()
            conn.close()

    def upload_property(self, struc_id: int, properties_dict: dict, project: str = "Property_Upload"):
        """
        Uploads new properties for an existing structure identified by struc_id.
        Creates a brand new upload_id in the Uploads table for this property upload event.
        
        properties_dict expects format: 
        {
            'property_name': {'value': float, 'unit': str, 'program': str}, ...
        }
        """
        if not properties_dict:
            print("No properties provided to upload.")
            return False

        username = os.environ.get("USER")
        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        try:
            # 1. Verify that the target structure exists
            cursor.execute("SELECT struc_id FROM Structures WHERE struc_id = %s;", (struc_id,))
            if not cursor.fetchone():
                print(f"Error: Structure with struc_id {struc_id} does not exist.")
                return False

            # 2. Create a new entry in the Uploads table
            upload_query = "INSERT INTO Uploads (username, project) VALUES (%s, %s);"
            cursor.execute(upload_query, (username, project))
            new_upload_id = cursor.lastrowid

            # 3. Insert the properties linked to struc_id and new_upload_id
            prop_query = """
            INSERT INTO Properties (struc_id, upload_id, property, value, unit, program)
            VALUES (%s, %s, %s, %s, %s, %s);
            """
            prop_payload = []
            for prop_name, data in properties_dict.items():
                prop_payload.append((
                    struc_id, 
                    new_upload_id, 
                    prop_name, 
                    data["value"], 
                    data.get("unit"), 
                    data.get("program")
                ))

            cursor.executemany(prop_query, prop_payload)

            # Commit the transaction
            conn.commit()
            print(f"Successfully uploaded {len(prop_payload)} propert(ies) for struc_id: {struc_id} under new upload_id: {new_upload_id}.")
            return new_upload_id

        except mysql.connector.Error as err:
            conn.rollback()
            print(f"Failed to upload properties, rolled back changes. Error: {err}")
            raise err
        finally:
            cursor.close()
            conn.close()

    def get_projects(self):
        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        query="""
            SELECT DISTINCT project FROM Uploads;
        """
        try:
            cursor.execute(query,)
            rows = cursor.fetchall()
            data=[]
            if not rows:
                return data
            
            for row in rows:
                data.append(row["project"])

        except mysql.connector.Error as err:
                    print(f"Database query failed. Error: {err}")
                    raise err
        finally:
            cursor.close()
            conn.close()
        return data
    
    def get_users(self):
        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        query="""
            SELECT DISTINCT username FROM Uploads;
        """
        try:
            cursor.execute(query,)
            rows = cursor.fetchall()
            data=[]
            if not rows:
                return data
            
            for row in rows:
                data.append(row["username"])

        except mysql.connector.Error as err:
                    print(f"Database query failed. Error: {err}")
                    raise err
        finally:
            cursor.close()
            conn.close()
        return data

    def fetch_table(self, table_name: str) -> list[dict]:
        """Fetch all rows from an allowed table."""
        ALLOWED_TABLES = {"Uploads", "Compositions", "Structures", "Properties"}
        if table_name not in ALLOWED_TABLES:
            raise ValueError(
                f"Invalid table name '{table_name}'. Allowed: {ALLOWED_TABLES}"
            )

        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        query = f"SELECT * FROM `{table_name}`;"

        try:
            cursor.execute(query)
            # fetchall() already returns a list of dicts (or empty list if no rows)
            return cursor.fetchall() or []

        except mysql.connector.Error as err:
            print(f"Database query failed. Error: {err}")
            raise err
        finally:
            cursor.close()
            conn.close()

    def get_property(self, by=None, entry=None):

        valid_criteria = ["comp", "composition", "project", "user"]
        if by not in valid_criteria or entry is None:
            print(f"Error: Invalid criteria. 'by' must be one of {valid_criteria} and 'entry' must be provided.")
            return False

        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        # SQL query to grab entry, structure info, and the specific energy property
        query = """
            SELECT 
                c.formula,
                s.struc_id,
                s.spacegroup,
                s.prototype,
                s.fu,
                s.natoms,
                s.dimension,
                p.property AS property_name,
                p.value,
                p.unit,
                p.program,
                p.prop_id,
                s.is_struct_new,
                s.is_SG_new,
                s.similar_to,
                s.relaxation,
                u.project,
                u.username,
                u.uploaded_at
            FROM Compositions c
            JOIN Structures s ON c.comp_id = s.comp_id
            JOIN Properties p ON s.struc_id = p.struc_id
            JOIN Uploads u ON s.upload_id = u.upload_id"""
        
        if((by=="comp") or (by=="composition")):
            entry=self.order_formula(entry)
            query+="\nWHERE c.formula = %s;"
        elif(by=="project"):
            query+="\nWHERE u.project = %s;"
        elif(by=="user"):
            query+="\nWHERE u.username = %s;"
        else:
            query+=";"
             
        try:
            if((by=="comp") or (by=="composition") or (by=="project") or (by=="user")):
                cursor.execute(query, (entry,))
            else:
                cursor.execute(query,)
            rows = cursor.fetchall()

            data=[]
            if not rows:
                return data

            for row in rows:
                # # 1. Decompress and deserialize the Pymatgen structure
                # decompressed_bytes = zlib.decompress(row['compressed_structure'])
                # structure_dict = json.loads(decompressed_bytes.decode('utf-8'))
                # struct_obj = Structure.from_dict(structure_dict)

                # 2. Build a data dictionary mimicking Materials Project format
                entry = {
                    "formula_pretty": row['formula'],
                    # "composition_reduced": struct_obj.composition.reduced_formula,
                    "struc_id": row['struc_id'],
                    "prop_id":row['prop_id'],
                    "property": row['property_name'],
                    "value": row['value'],
                    "unit": row['unit'],
                    "program": row['program'],
                    "spacegroup": row['spacegroup'],
                    "prototype": row['prototype'],
                    "fu": row['fu'],
                    "natoms": row['natoms'],
                    "dimension": row['dimension'],
                    "is_struct_new": row['is_struct_new'],
                    "is_SG_new": row['is_SG_new'],
                    "similar_to": row['similar_to'],
                    "relaxation": row['relaxation'],
                    "project": row['project'],
                    "user": row['username'],
                    "upload_time": row['uploaded_at']
                }
                data.append(entry)
            return data

        except mysql.connector.Error as err:
                    print(f"Database query failed. Error: {err}")
                    raise err
        finally:
            cursor.close()
            conn.close()


    def get_structure(self, by=None, entry=None):

        valid_criteria = ["comp", "composition", "project", "user"]
        if by not in valid_criteria or entry is None:
            print(f"Error: Invalid criteria. 'by' must be one of {valid_criteria} and 'entry' must be provided.")
            return False

        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        # SQL query to grab entry, structure info, and the specific energy property
        query = """
            SELECT 
                c.formula,
                s.struc_id,
                s.structure AS compressed_structure,
                s.spacegroup,
                s.prototype,
                s.fu,
                s.natoms,
                s.dimension,
                s.is_struct_new,
                s.is_SG_new,
                s.similar_to,
                s.relaxation,
                u.project,
                u.username,
                u.uploaded_at
            FROM Compositions c
            JOIN Structures s ON c.comp_id = s.comp_id
            JOIN Uploads u ON s.upload_id = u.upload_id"""
        
        if((by=="comp") or (by=="composition")):
            entry=self.order_formula(entry)
            query+="\nWHERE c.formula = %s;"
        elif(by=="project"):
            query+="\nWHERE u.project = %s;"
        elif(by=="user"):
            query+="\nWHERE u.username = %s;"
        else:
            query+=";"
             
        try:
            if((by=="comp") or (by=="composition") or (by=="project") or (by=="user")):
                cursor.execute(query, (entry,))
            else:
                cursor.execute(query,)
            rows = cursor.fetchall()

            data=[]
            if not rows:
                return data

            for row in rows:
                # 1. Decompress and deserialize the Pymatgen structure
                decompressed_bytes = zlib.decompress(row['compressed_structure'])
                structure_dict = json.loads(decompressed_bytes.decode('utf-8'))
                struct_obj = Structure.from_dict(structure_dict)

                # 2. Build a data dictionary mimicking Materials Project format
                entry = {
                    "formula_pretty": row['formula'],
                    "struc_id": row['struc_id'],
                    "structure": struct_obj,
                    "spacegroup": row['spacegroup'],
                    "prototype": row['prototype'],
                    "fu": row['fu'],
                    "natoms": row['natoms'],
                    "dimension": row['dimension'],
                    "is_struct_new": row['is_struct_new'],
                    "is_SG_new": row['is_SG_new'],
                    "similar_to": row['similar_to'],
                    "relaxation": row['relaxation'],
                    "project": row['project'],
                    "user": row['username'],
                    "upload_time": row['uploaded_at']
                }
                data.append(entry)
            return data

        except mysql.connector.Error as err:
                    print(f"Database query failed. Error: {err}")
                    raise err
        finally:
            cursor.close()
            conn.close()
    
    def write_cif(self,struc_id,path):
        from pymatgen.io.cif import CifWriter
        import os

        conn = self._get_connection()
        cursor = conn.cursor(dictionary=True)

        # SQL query to grab entry, structure info, and the specific energy property
        query = """
            SELECT 
            c.formula,
            s.spacegroup,
            s.structure AS compressed_structure
            FROM Structures s
            JOIN Compositions c ON c.comp_id = s.comp_id
            WHERE struc_id = %s;"""
        
        try:
            cursor.execute(query, (struc_id,))
            row = cursor.fetchall()[0]
            if not row:
                return row
            
            os.makedirs(path, exist_ok=True)
            formula= row['formula']
            spacegroup= row['spacegroup']
            decompressed_bytes = zlib.decompress(row['compressed_structure'])
            structure_dict = json.loads(decompressed_bytes.decode('utf-8'))
            struct_obj = Structure.from_dict(structure_dict)

            struc_path=os.path.join(path,formula+"_sym"+str(spacegroup)+"_s"+str(struc_id)+".cif")
            writer = CifWriter(struct_obj, symprec=1e-2, angle_tolerance=5, refine_struct=True)
            writer.write_file(struc_path)
            return struc_path

        except mysql.connector.Error as err:
                    print(f"Database query failed. Error: {err}")
                    raise err
        finally:
            cursor.close()
            conn.close()

    def generate_castep_input(self,struc_id,path):
        import autocasp

        struc_path=self.write_cif(struc_id,path)
        autocasp.run(cif_file=struc_path)
    
    def clear_database(self):
            """
            Safely clears all data from Uploads, Structures, Properties, and Compositions,
            resetting all auto-increment IDs back to 1.
            """

            choice = input("Do you want to delete all entries in the database? (y/n): ").strip().lower()
            if choice!="y":
                 print("Operation aborted.")
                 return False
            
            conn = self._get_connection()
            cursor = conn.cursor()

            try:
                print("Clearing all data from TCCDEM_DB...")
                
                # 1. Disable foreign key checks temporarily to allow TRUNCATE
                cursor.execute("SET FOREIGN_KEY_CHECKS = 0;")
                
                # 2. Clear tables and reset auto-increment counters
                tables_to_clear = ["Properties", "Structures", "Compositions", "Uploads"]
                for table in tables_to_clear:
                    cursor.execute(f"TRUNCATE TABLE {table};")
                    print(f" - Emptied table: {table}")
                    
                # 3. Re-enable foreign key checks
                cursor.execute("SET FOREIGN_KEY_CHECKS = 1;")
                
                conn.commit()
                print("Database successfully cleared and reset.")
                
            except mysql.connector.Error as err:
                conn.rollback()
                print(f"Failed to clear database. Error: {err}")
                raise err
            finally:
                cursor.close()
                conn.close()

    def clear_entries(self, by=None, entry=None):
        """
        Deletes specific entries from database tables based on criteria:
        Parameters:
            by (str): 'struc_id', 'prop_id', 'project', or 'user' (or 'username')
            entry (str | int): The value matching the 'by' filter.
        """
        valid_criteria = ["struc_id", "prop_id", "project", "user", "username"]
        if by not in valid_criteria or entry is None:
            print(f"Error: Invalid criteria. 'by' must be one of {valid_criteria} and 'entry' must be provided.")
            return False

        # Normalize username/user key
        filter_by = "username" if by == "user" else by

        choice = input(f"Are you sure you want to delete entries where {filter_by} = '{entry}'? (y/n): ").strip().lower()
        if choice != "y":
            print("Operation aborted.")
            return False

        conn = self._get_connection()
        cursor = conn.cursor()

        try:
            # Case 1: Delete only a single Property record
            if filter_by == "prop_id":
                cursor.execute("DELETE FROM Properties WHERE prop_id = %s;", (entry,))
                deleted_rows = cursor.rowcount
                print(f"Deleted {deleted_rows} record from Properties.")

            # Case 2: Delete by Structure ID (Deletes child Properties, then Structure)
            elif filter_by == "struc_id":
                # Delete associated properties first
                cursor.execute("DELETE FROM Properties WHERE struc_id = %s;", (entry,))
                prop_count = cursor.rowcount

                # Delete the structure itself
                cursor.execute("DELETE FROM Structures WHERE struc_id = %s;", (entry,))
                struc_count = cursor.rowcount
                print(f"Deleted {struc_count} structure and {prop_count} associated property record(s).")

            # Case 3: Delete by Project or User
            elif filter_by in ["project", "username"]:
                column_name = "project" if filter_by == "project" else "username"
                
                # Fetch target upload IDs
                cursor.execute(f"SELECT upload_id FROM Uploads WHERE {column_name} = %s;", (entry,))
                upload_rows = cursor.fetchall()
                
                if not upload_rows:
                    print(f"No records found for {filter_by} = '{entry}'.")
                    return False
                
                upload_ids = [row[0] for row in upload_rows]
                format_strings = ','.join(['%s'] * len(upload_ids))

                # 1. Delete associated Properties
                cursor.execute(f"DELETE FROM Properties WHERE upload_id IN ({format_strings});", tuple(upload_ids))
                prop_count = cursor.rowcount

                # 2. Delete associated Structures
                cursor.execute(f"DELETE FROM Structures WHERE upload_id IN ({format_strings});", tuple(upload_ids))
                struc_count = cursor.rowcount

                # 3. Delete Upload records
                cursor.execute(f"DELETE FROM Uploads WHERE upload_id IN ({format_strings});", tuple(upload_ids))
                upload_count = cursor.rowcount

                print(f"Deleted {upload_count} upload(s), {struc_count} structure(s), and {prop_count} property record(s).")

            # Clean up orphaned compositions that no longer have any associated structures
            cursor.execute("""
                DELETE FROM Compositions 
                WHERE comp_id NOT IN (SELECT DISTINCT comp_id FROM Structures WHERE comp_id IS NOT NULL);
            """)
            orphaned_comps = cursor.rowcount # returns the number of rows affected or retrieved by the most recently executed SQL statement
            if orphaned_comps > 0:
                print(f"Cleaned up {orphaned_comps} orphaned composition record(s).")

            conn.commit()
            print("Deletion completed successfully.")
            return True

        except mysql.connector.Error as err:
            conn.rollback()
            print(f"Failed to delete entries. Changes rolled back. Error: {err}")
            raise err
        finally:
            cursor.close()
            conn.close()