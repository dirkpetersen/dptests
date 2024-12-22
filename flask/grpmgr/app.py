from flask import Flask, render_template, jsonify, request
from datetime import datetime
from ldap3 import Server, Connection, SUBTREE
from dotenv import load_dotenv
import os
import json

# Load environment variables
load_dotenv()

app = Flask(__name__)

def get_ldap_connection():
    host = os.getenv('LDAP_HOST', 'localhost')
    port = int(os.getenv('LDAP_PORT', '10389'))
    bind_dn = os.getenv('LDAP_BIND_DN', 'uid=svc-dri-lookup,ou=Specials,o=orst.edu')
    password = os.getenv('LDAP_PASSWORD')
    
    if not password:
        raise ValueError("LDAP_PASSWORD environment variable must be set")
        
    server = Server(host=host, port=port)
    return Connection(
        server,
        user=bind_dn,
        password=password,
        auto_bind=True
    )

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/users/search')
def search_users():
    query = request.args.get('q', '')
    if not query:
        return jsonify([])
    
    ldap_filter = f"(&(objectClass=person)(|(cn=*{query}*)(mail=*{query}*)(ou=*{query}*)))"
    
    with get_ldap_connection() as conn:
        conn.search(
            os.getenv('LDAP_BASE_DN'),
            ldap_filter,
            SUBTREE,
            attributes=['uid', 'cn', 'mail', 'ou']
        )
        
        results = []
        for entry in conn.entries:
            results.append({
                "id": entry.uid.value,
                "name": entry.cn.value,
                "email": entry.mail.value if hasattr(entry, 'mail') else '',
                "department": entry.ou.value if hasattr(entry, 'ou') else ''
            })
            
    return jsonify(results)

@app.route('/api/groups/search')
def search_groups():
    query = request.args.get('q', '')
    if not query:
        return jsonify([])
    
    ldap_filter = f"(&(objectClass=groupOfNames)(cn=*{query}*)))"
    
    with get_ldap_connection() as conn:
        conn.search(
            os.getenv('LDAP_BASE_DN'),
            ldap_filter,
            SUBTREE,
            attributes=['cn', 'description']
        )
        
        results = []
        for entry in conn.entries:
            results.append({
                "id": entry.cn.value,
                "name": entry.cn.value,
                "description": entry.description.value if hasattr(entry, 'description') else ''
            })
            
    return jsonify(results)

@app.route('/api/submit-changes', methods=['POST'])
def submit_changes():
    data = request.json
    group_id = data.get('group')
    user_ids = data.get('users', [])
    
    if not group_id or not user_ids:
        return jsonify({
            "status": "error",
            "message": "Missing group or users"
        }), 400
    
    try:
        with get_ldap_connection() as conn:
            # Find the group DN
            group_filter = f"(&(objectClass=groupOfNames)(cn={group_id}))"
            conn.search(
                os.getenv('LDAP_BASE_DN'),
                group_filter,
                SUBTREE,
                attributes=['member']
            )
            
            if not conn.entries:
                return jsonify({
                    "status": "error",
                    "message": "Group not found"
                }), 404
                
            group_dn = conn.entries[0].entry_dn
            
            # Add users to group
            for user_id in user_ids:
                user_dn = f"uid={user_id},{os.getenv('LDAP_BASE_DN')}"
                conn.modify(group_dn, {'member': [(2, [user_dn])]})  # 2 = MODIFY_ADD
                
        return jsonify({
            "status": "success",
            "message": "Changes submitted successfully",
            "timestamp": datetime.now().isoformat()
        })
        
    except Exception as e:
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5555,debug=True)
