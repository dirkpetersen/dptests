from flask import Flask, render_template, jsonify, request
from datetime import datetime
import json

app = Flask(__name__)

# Mock data - in production this would come from a database
USERS = [
    {"id": 1, "name": "John Doe", "email": "john@example.com", "department": "IT"},
    {"id": 2, "name": "Jane Smith", "email": "jane@example.com", "department": "HR"},
    {"id": 3, "name": "Bob Wilson", "email": "bob@example.com", "department": "Finance"},
    {"id": 4, "name": "Alice Brown", "email": "alice@example.com", "department": "IT"},
]

GROUPS = [
    {"id": 1, "name": "IT Admin", "description": "IT Administrators"},
    {"id": 2, "name": "HR Team", "description": "Human Resources Team"},
    {"id": 3, "name": "Finance Users", "description": "Finance Department Users"},
]

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/users/search')
def search_users():
    query = request.args.get('q', '').lower()
    results = [user for user in USERS if query in user['name'].lower() or 
               query in user['email'].lower() or 
               query in user['department'].lower()]
    return jsonify(results)

@app.route('/api/groups/search')
def search_groups():
    query = request.args.get('q', '').lower()
    results = [group for group in GROUPS if query in group['name'].lower() or 
               query in group['description'].lower()]
    return jsonify(results)

@app.route('/api/submit-changes', methods=['POST'])
def submit_changes():
    data = request.json
    # In a real application, you would process the changes here
    # For now, we'll just return a success message
    return jsonify({
        "status": "success",
        "message": "Changes submitted successfully",
        "timestamp": datetime.now().isoformat()
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5555,debug=True)
