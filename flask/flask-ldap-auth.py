#!/usr/bin/env python3

from flask import Flask, Response, request, render_template_string
import base64, ldap, logging

app = Flask(__name__)
app.logger.setLevel(logging.INFO)

LDAPHOST = 'ds2.acc.xxxx.xxx'
UOUSERS = 'ou=Users,dc=acc,dc=xxxx,dc=xxx'
ALLOWED_USERS = ['bedricks', 'peterdir']  # List of allowed usernames

login_html = '''
<!DOCTYPE html>
<html>
<head>
    <title>Login</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            background-color: #f2f2f2;
        }
        .login-form {
            width: 300px;
            margin: 0 auto;
            padding: 20px;
            background-color: #fff;
            border-radius: 5px;
            box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
        }
        .login-form h2 {
            text-align: center;
            margin-bottom: 20px;
        }
        .login-form input[type="text"],
        .login-form input[type="password"] {
            width: 100%;
            padding: 10px;
            margin-bottom: 10px;
            border: 1px solid #ccc;
            border-radius: 3px;
        }
        .login-form input[type="submit"] {
            width: 100%;
            padding: 10px;
            background-color: #4CAF50;
            color: #fff;
            border: none;
            border-radius: 3px;
            cursor: pointer;
        }
    </style>
</head>
<body>
    <div class="login-form">
        <h2>Please enter your OHSU username and password:</h2>
        <form method="post" autocomplete="off">
            <input type="text" name="username" placeholder="Username" required autocomplete="off">
            <input type="password" name="password" placeholder="Password" required autocomplete="off">
            <input type="submit" value="Login">
        </form>
    </div>
</body>
</html>
'''

@app.route("/", methods=['GET', 'POST'])
def default():
    msg = 'Wrong URL! You must access /pyauth to authenticate !'
    app.logger.info(msg)
    return Response(msg, 404)

@app.route("/pyauth", methods=['GET', 'POST'])
def auth():
    auth_header = request.headers.get('Authorization')
    if not auth_header:
        if request.method == 'POST':
            username = request.form['username']
            password = request.form['password']
            if ldap_authenticate(username, password):
                if username in ALLOWED_USERS:
                    app.logger.info(f"Authenticated user: {username}")
                    return Response('Authenticated', 200)
                else:
                    app.logger.info(f"Authorization failed for user: {username}")
                    return Response('You are not a member of allowed_users', 403)
            else:
                app.logger.info(f"Authentication failed for user: {username}")
                return Response('Unauthorized', 401)
        return render_template_string(login_html)

    username, password = parse_auth_header(auth_header)
    if username not in ALLOWED_USERS:
        app.logger.info(f"Authorization failed: {username} is not in ALLOWED_USERS")
        return Response('You are not a member of ALLOWED_USERS', 403)  # 403 Forbidden Access

    if ldap_authenticate(username, password):
        app.logger.info(f"Authenticated user: {username}")
        return Response('Authenticated', 200)
    else:
        app.logger.info(f"Authentication failed for user: {username}")
        return Response('Unauthorized', 401)

def parse_auth_header(header):
    header = header.replace("Basic ", "", 1)
    username, password = base64.b64decode(header).decode('utf-8').split(':')
    return username, password

def ldap_authenticate(username, password):
    try:
        conn = ldap.initialize(f'ldap://{LDAPHOST}')  # ldaps is not working
        conn.simple_bind_s(f'uid={username},{UOUSERS}', password)
        return True
    except ldap.LDAPError:
        return False

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8999)
