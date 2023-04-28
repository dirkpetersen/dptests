import os, sys
from flask import Flask, redirect, url_for, render_template, flash
from flask_dance.contrib.google import make_google_blueprint, google
from flask_dance.contrib.github import make_github_blueprint, github
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager, UserMixin, login_user, logout_user, login_required, current_user

app = Flask(__name__)
app.config["SECRET_KEY"] = os.environ.get("FLASK_SECRET_KEY")
app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///users.db"
db = SQLAlchemy(app)

login_manager = LoginManager(app)

class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(80), unique=True, nullable=False)
    provider = db.Column(db.String(10), nullable=False)

@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

google_blueprint = make_google_blueprint(
    client_id=os.environ.get("GOOGLE_OAUTH_CLIENT_ID"),
    client_secret=os.environ.get("GOOGLE_OAUTH_CLIENT_SECRET"),
    scope=["profile", "email"]
)

github_blueprint = make_github_blueprint(
    client_id="your_github_client_id",
    client_secret="your_github_client_secret",
    scope="user:email"
)

app.register_blueprint(google_blueprint, url_prefix="/login/google")
app.register_blueprint(github_blueprint, url_prefix="/login/github")

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/login")
def login():
    return render_template("login.html")

@app.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for("index"))

@app.route("/google_authorized")
def google_authorized():
    if not google.authorized:
        return redirect(url_for("google.login"))

    resp = google.get("/oauth2/v2/userinfo")
    if resp.ok:
        user_data = resp.json()
        user = User.query.filter_by(username=user_data["email"], provider="google").first()
        if not user:
            user = User(username=user_data["email"], provider="google")
            db.session.add(user)
            db.session.commit()
        login_user(user)
        return redirect(url_for("index"))

    flash("Failed to fetch user data", "error")
    return redirect(url_for("index"))

@app.route("/github_authorized")
def github_authorized():
    if not github.authorized:
        return redirect(url_for("github.login"))

    resp = github.get("/user")
    if resp.ok:
        user_data = resp.json()
        user = User.query.filter_by(username=user_data["login"], provider="github").first()
        if not user:
            user = User(username=user_data["login"], provider="github")
            db.session.add(user)
            db.session.commit()
        login_user(user)
        return redirect(url_for("index"))

    flash("Failed to fetch user data", "error")
    return redirect(url_for("index"))

if __name__ == "__main__":
    app.run(debug=True)
