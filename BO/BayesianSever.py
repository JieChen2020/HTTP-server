import requests
import json

# API的基本URL
base_url = 'http://10.97.40.181:123'

# 初始请求数据
data = {
    "condition": [[11, 6, 'Pd(PPh3)4', 'Na2CO3'], [4.00, 5.72, 'SPhos Pd G3', 'K3PO4'], [7.48, 2.13, 'PdCl2(dppf)', 'K2CO3'],
                  [6.21, 5.21, 'PdCl2(dppf)', 'Na2CO3'], [2.60, 4.08, 'XPhos Pd G3', 'K2CO3']],
    "outcomes": [0.164, 0.127, 0.538, 0.541, 0.101],
    "q": 4,
    "goal": [
        {
            "name": "yield",
            "target": "max"
        }
    ],
    "design_space": [
        {
            "name": "Catalyst loading (mol %)",
            "type": "continuous",
            "range": [
                2,
                12
            ]
        },
        {
            "name": "Reaction time (hr)",
            "type": "continuous",
            "range": [
                2,
                6
            ]
        },
        {
            "name": "Catalyst system",
            "type": "categorical",
            "range": [
                'Pd(PPh3)4',
                'PdCl2(dppf)',
                'XPhos Pd G3',
                'SPhos Pd G3'
            ]
        },
        {
            "name": "Base",
            "type": "categorical",
            "range": [
                'Na2CO3',
                'K2CO3',
                'K3PO4'
            ]
        }
    ],
    "num_of_init": 4,
    "reaction": "aa"
}

# 将字典转换为JSON格式
json_data = json.dumps(data)

# 设置请求的头部，表明我们发送的是JSON格式
headers = {'Content-Type': 'application/json'}

# 发送POST请求到/get-next-exeps
response_get_next_exeps = requests.post(f'{base_url}/get-next-exps', headers=headers, data=json_data)
print("Status code:", response_get_next_exeps.status_code)
print("Response text:", response_get_next_exeps.text)
if response_get_next_exeps.status_code == 200:
    print("Response from /get-next-exps:")
    print(response_get_next_exeps.json())
else:
    print("Server error!")
