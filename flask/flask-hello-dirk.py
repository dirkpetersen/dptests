#! /usr/bin/env python3

from flask import Flask, render_template_string

app = Flask(__name__)

@app.route('/')
def hello():
    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Hello Dirk</title>
        <style>
            body {
                font-family: Arial, sans-serif;
                text-align: center;
                background-color: #f0f0f0;
            }
            h1 {
                color: #333;
                font-size: 3em;
                margin-top: 50px;
            }
            p {
                color: #666;
                font-size: 1.5em;
            }
        </style>
    </head>
    <body>
        <h1>Hello Dirk!</h1>
        <p>Welcome to this simple Flask app.</p>
    </body>
    </html>
    """
    return render_template_string(html)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8998)