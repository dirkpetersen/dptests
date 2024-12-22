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
    ldap_url = os.getenv('LDAP_URL', 'ldap://localhost:10389')
    bind_dn = os.getenv('LDAP_BIND_DN', 'uid=svc-dri-lookup,ou=Specials,o=orst.edu')
    password = os.getenv('LDAP_PASSWORD')
    
    if not password:
        raise ValueError("LDAP_PASSWORD environment variable must be set")
        
    server = Server(ldap_url)
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
    query = request.args.get('q', '').strip()
    if not query:
        return jsonify([])
    
    # Split the query into words
    words = query.split()
    
    # Create an AND filter that requires all words to be present in any field
    word_filters = []
    for word in words:
        # Check if word is numeric for uidNumber search
        is_numeric = word.isdigit()
        # Include uidNumber only if the word is numeric
        base_filter = f"(|(cn=*{word}*)(mail=*{word}*)(uid=*{word}*)(title=*{word}*)(osuDepartment=*{word}*))"
        word_filter = f"(|{base_filter}(uidNumber={word}))" if is_numeric else base_filter
        word_filters.append(word_filter)
    
    # Combine with AND to require all words
    ldap_filter = f"(&(objectClass=person){''.join(word_filters)})"
    
    with get_ldap_connection() as conn:
        conn.search(
            os.getenv('LDAP_BASE_DN_USER'),
            ldap_filter,
            SUBTREE,
            attributes=['uid', 'cn', 'mail', 'title', 'osuDepartment', 'uidNumber']
        )
        
        results = []
        for entry in conn.entries:
            if (hasattr(entry, 'mail') and entry.mail.value and 
                hasattr(entry, 'uid') and entry.uid.value and 
                entry.uid.value.lower() not in ['null', 'none']):  # Only include entries with email and valid UID
                results.append({
                    "id": entry.uid.value,
                    "name": entry.cn.value,
                    "email": entry.mail.value,
                    "title": entry.title.value if hasattr(entry, 'title') and entry.title.value not in ['null', 'NULL', None] else '',
                    "department": entry.osuDepartment.value if hasattr(entry, 'osuDepartment') and entry.osuDepartment.value not in ['null', 'NULL', None] else ''
                })
            
    return jsonify(results)

@app.route('/api/groups/search')
def search_groups():
    query = request.args.get('q', '').strip()
    if not query:
        return jsonify([])
    
    # Split the query into words
    words = query.split()
    word_filters = []
    for word in words:
        # Check if word is numeric for gidNumber search
        is_numeric = word.isdigit()
        # Include gidNumber only if the word is numeric
        if is_numeric:
            word_filter = f"(gidNumber={word})"
        else:
            word_filter = f"(cn=*{word}*)"
        word_filters.append(word_filter)
    
    # Combine with AND to require all words
    ldap_filter = f"(&(objectClass=posixGroup){''.join(word_filters)})"
    
    with get_ldap_connection() as conn:
        conn.search(
            "ou=group,o=orst.edu",
            ldap_filter,
            SUBTREE,
            attributes=['cn', 'gidNumber', 'memberUid']
        )
        
        results = []
        for entry in conn.entries:
            results.append({
                "id": entry.cn.value,
                "name": entry.cn.value,
                "description": "",
                "gidNumber": entry.gidNumber.value if hasattr(entry, 'gidNumber') else None,
                "members": entry.memberUid.values if hasattr(entry, 'memberUid') else []
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
            group_filter = f"(&(objectClass=posixGroup)(cn={group_id}))"
            conn.search(
                os.getenv('LDAP_BASE_DN_GROUP'),
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
                user_dn = f"uid={user_id},{os.getenv('LDAP_BASE_DN_USER')}"
                conn.modify(group_dn, {'memberUid': [(2, [user_id])]})  # 2 = MODIFY_ADD
                
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
