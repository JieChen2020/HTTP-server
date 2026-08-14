from flask import Flask, request, jsonify
import requests
import json
import pandas as pd
import os
from anytree import Node, RenderTree
import copy
from rdkit import Chem
from rdkit.Chem import Descriptors


app = Flask(__name__)
EXCEL_PATH = 'smiles_taskId.xlsx'
EXCEL_PATH_2 = 'smiles_price.xlsx'
URL = "https://ai-cn.chemlex.com/retroplanning"
HEADERS = {
    "X-SECRET-API-KEY": "96667896156d976a3ab95d2a42e6a3f302221957e8600f1e342864b593d9f4fcaa756fb13beab1dc549bb88d3340e786526f8f8eed666923077a0389d7c495c0",
    "Content-Type": "application/json",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
}
current_task_id = 1


@app.route('/ping', methods=['GET'])
def ping():
    # Http server status
    return jsonify({"status": "online"})


@app.route('/receive_data', methods=['POST'])
def receive_data():
    global current_task_id

    try:
        request_data = request.get_json()
        data = request_data.get('data')
        print(data)

        if data != "query":
            tool_input_str = json.loads(data).get('SMILES')
            print(tool_input_str)
            df = load_excel()
            task_row = df[df['smiles'] == tool_input_str]
            if not task_row.empty:
                current_task_id = int(task_row.iloc[0]['taskId'])
                print(current_task_id)
                response = ({"status": "running"})
            else:
                current_task_id = submit_task(tool_input_str)
                if isinstance(current_task_id, int):
                    save_task_to_excel(tool_input_str, current_task_id)
                    response = {"status": "running"}
                else:
                    current_task_id = 1  # initialize
                    response = {"status": "error", "message": "Submit task failed"}

        else:
            print(current_task_id)
            if current_task_id == 1:
                response = {"status": "error", "message": "No valid taskId available"}
            else:
                result = get_result(current_task_id)
                if result.get("status") == "completed":
                    result["status"] = "success"
                    response = result
                    print(response)
                    current_task_id = 1
                else:
                    response = result

        return jsonify(response)

    except Exception as e:
        current_task_id = 1
        return jsonify({"status": "error", "message": str(e)})


def submit_task(smiles):
    try:
        data = {
            "method": "submit_task",
            "param": {
                "param":
                    {
                        "smiles": smiles,
                        "param": {
                            "search_mode": {
                                "use_expert_rules": False,
                                "use_model": True,
                                "creative_mode": False,
                                "fast_mode": True
                            },
                            "max_depth": 3,
                            "disliked_reactions": [],
                            "disliked_compounds": [],
                            "max_simulation_time": 60
                        }
                    }
            }
        }
        payload = json.dumps(data)
        response = requests.post(URL, headers=HEADERS, data=payload)
        if response.status_code == 200:
            return response.json().get("data")
        else:
            return None
    except Exception as e:
        return None


def get_result(taskId):
    try:
        data = {
            "method": "get_result",
            "param": {"taskId": taskId}
        }
        payload = json.dumps(data)
        response = requests.post(URL, headers=HEADERS, data=payload)

        if response.status_code == 200:
            res_json = json.loads(response.text)
            status = res_json["data"].get("status", "")
            process = res_json["data"].get("process", 0)

            if status == "completed":
                result_raw = res_json["data"].get("result")
                if not result_raw:
                    return {"status": "error", "message": "Result field is empty"}
                try:
                    list_of_paths = json.loads(result_raw)
                    result_dict = {"status": "completed"}
                    j = 1
                    for i in range(50):
                        if list_of_paths[i]['is_buyable_route'] and j < 50:
                            path_info = list_of_paths[i]['path']
                            print(path_info)
                            simplified_path_info = extract_relevant_fields(path_info)
                            result_dict[f"path{j}_info"] = copy.deepcopy(simplified_path_info)
                            reaction_tree = extract_reaction_tree_dict(simplified_path_info)
                            result_dict[f"path{j}"] = reaction_tree
                            j = j + 1

                    return result_dict
                except Exception as e:
                    return {"status": "error", "message": f"Failed to parse result: {str(e)}"}

            elif status in ["pending", "start"]:
                return {"status": "running", "message": f"Task is {status}. Current process: {process}%"}

            elif status == "fail":
                return {"status": "failed", "message": "Task failed during execution"}

            else:
                return {"status": "error", "message": f"Unknown status: {status}"}

        else:
            return {"status": "error", "message": f"HTTP {response.status_code}"}
    except Exception as e:
        return f"Error: {str(e)}"


def save_task_to_excel(smiles, taskId):
    df = load_excel()
    new_row = pd.DataFrame([[smiles, taskId]], columns=['smiles', 'taskId'])
    df = pd.concat([df, new_row], ignore_index=True)
    df.to_excel(EXCEL_PATH, index=False)


def load_excel():
    if os.path.exists(EXCEL_PATH):
        return pd.read_excel(EXCEL_PATH)
    else:
        return pd.DataFrame(columns=['smiles', 'taskId'])


def load_excel_price():
    if os.path.exists(EXCEL_PATH_2):
        return pd.read_excel(EXCEL_PATH_2)
    else:
        return pd.DataFrame(columns=['smiles', 'price'])


def extract_relevant_fields(node: dict) -> dict:
    """
    Recursively extracts only the required fields from a reaction tree node:
    - smiles
    - terminal
    - children (recursively)
    - id
    - price
    - reaction_conditions (effective)
    """
    def extract_valid_condition(cond: dict) -> dict:
        def flatten_entries(source):
            """Flatten source data into a list of entries, compatible with both [[entry, ...]] and [entry, ...] structures"""
            entries = []
            if not source:
                return entries
            for item in source:
                if isinstance(item, list):
                    entries.extend(item)
                elif isinstance(item, dict):
                    entries.append(item)
            return entries

        # Merge from_patents and from_expert
        patents = cond.get("from_patents", [])
        experts_raw = cond.get("from_experts", [])
        experts = [e.get("expert_recommended", e) if isinstance(e, dict) else e for e in (experts_raw if isinstance(experts_raw, list) else [experts_raw])]
        all_entries = flatten_entries(patents) + flatten_entries(experts)
        if all_entries:
            best_entry = {}
            max_yield = 30  # default yield is 30 if "yield" is "N/A"

            for entry in all_entries:
                preparation = entry.get("preparation", "N/A")
                if preparation == "N/A":
                    continue

                # Default "yield" to 30 if it's "N/A"
                entry_yield = entry.get("yield", "N/A")
                entry_yield = 30 if entry_yield == "N/A" else float(entry_yield)

                # Check if this entry has a higher yield
                if entry_yield > max_yield:
                    max_yield = entry_yield
                    best_entry = entry

            if best_entry:
                return best_entry

        return {}

    result = {}

    if 'smiles' in node:
        result['smiles'] = node['smiles']
    if 'id' in node:
        result['id'] = node['id']
    if 'terminal' in node:
        result['terminal'] = node['terminal']
    if 'compound_attributes' in node:
        result['price'] = node['compound_attributes']['price']

    if 'reaction_conditions' in node:
        selected_condition = extract_valid_condition(node['reaction_conditions'])
        result['reaction_conditions'] = selected_condition

    if 'children' in node:
        result['children'] = [extract_relevant_fields(child) for child in node['children']]
    else:
        result['children'] = []

    return result


def extract_reaction_tree_dict(data_dict: dict) -> dict:
    df = load_excel_price()
    level_info = {}

    def get_price(smiles_str: str):
        row = df[df['smiles'] == smiles_str]
        if not row.empty and 'price' in row.columns:
            return row.iloc[0]['price']
        else:
            return 'N.A.'

    def get_name(smiles_str: str):
        row = df[df['smiles'] == smiles_str]
        if not row.empty and 'name' in row.columns:
            return row.iloc[0]['name']
        else:
            return 'N.A.'

    def get_mw(smiles_str: str):
        """Calculate molecular weight from SMILES using RDKit, return 'N.A.' on failure"""
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is None:
                return 'N.A.'
            return round(Descriptors.MolWt(mol), 1)
        except Exception:
            return 'N.A.'

    def get_reagent_info(smiles_str: str):
        """Get name, price, and molecular weight simultaneously, return [name, price, mw]"""
        price = get_price(smiles_str)
        name = get_name(smiles_str)
        mw = get_mw(smiles_str)
        return {
                "Reagent SMILES": smiles_str,
                "Name": name,
                "Price": price,
                "Molecular Weight": mw,
            }

    def traverse(node, level=0):
        if node.get("terminal", True):
            return

        children = node.get("children", []) or []

        reactant = []
        for child in children:
            smi = child.get("smiles", "")
            reactant.append(get_reagent_info(smi))

        product_smiles = node.get("smiles", "")
        product = get_reagent_info(product_smiles)

        # Only record the first time this level is encountered
        if level not in level_info:
            reaction_conditions = node.get("reaction_conditions", {}) or {}
            reaction_conditions = dict(reaction_conditions)  # shallow copy

            level_info[level] = {
                "Product": product,
                "Reactant": reactant,
                "Synthetic description for reference": reaction_conditions.get('preparation', 'N.A.'),
            }

        # Recurse into child nodes
        for child in children:
            traverse(child, level + 1)

        # Traverse the entire tree first

    traverse(data_dict)

    result = {}

    if level_info:
        # Maximum level, e.g. Level 0,1,2 -> max_level = 2
        max_level = max(level_info.keys())

        # Map level to Step: Level 0 -> Step N (max), Level 1 -> Step N-1, ..., Level N -> Step 1
        for level in sorted(level_info.keys()):
            step_index = max_level - level + 1
            step_key = f"Step {step_index}"
            result[step_key] = level_info[level]

    return result


if __name__ == '__main__':
    app.run(host='', port=8002)
    # print(get_result(21))
    # print(submit_task("FC(F)(F)C1=CC(=NN1C1=CC=C(Br)C=C1)C1=CC=CC=C1"))
