#! /usr/bin/env python3

import os, sys
from flask import Flask, render_template, redirect, url_for, flash, request
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
#from flask_dance.contrib.google import make_google_blueprint, google
#from flask_dance.consumer.storage.sqla import OAuthConsumerMixin, SQLAlchemyStorage
#from flask_dance.consumer.backend.sqla import OAuthConsumerMixin, SQLAlchemyBackend
from flask_dance.contrib.google import make_google_blueprint, google
#from flask_dance.contrib.github import make_github_blueprint, github
#from flask_dance.consumer.backend.sqla import OAuthConsumerMixin, SQLAlchemyBackend
from flask_dance.consumer.storage.sqla import SQLAlchemyStorage
from flask_dance.consumer import oauth_authorized, oauth_error


from sqlalchemy.orm.exc import NoResultFound

# Requirements: 
# python3 -m pip install flask flask_sqlalchemy flask_login flask_dance[sqla] flask-dance[google]

app = Flask(__name__)

# Configurations
app.config['SECRET_KEY'] = os.environ.get("FLASK_SECRET_KEY")
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///q_and_a.sqlite'

# Models
class User(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100))
    email = db.Column(db.String(100), unique=True)
    questions = db.relationship('Question', backref='user', lazy=True)

class Question(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    data = db.Column(db.Text)
    response = db.Column(db.Text, nullable=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))

class OAuth(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    provider = db.Column(db.String(50), nullable=False)
    provider_user_id = db.Column(db.String(256), nullable=False)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'), nullable=False)
    user = db.relationship(User)

# Initialize
db = SQLAlchemy(app)
login_manager = LoginManager(app)

# Google OAuth Blueprint
blueprint = make_google_blueprint(
    client_id=os.environ.get("GOOGLE_OAUTH_CLIENT_ID"),
    client_secret=os.environ.get("GOOGLE_OAUTH_CLIENT_SECRET"),
    scope=["profile", "email"],
    redirect_url="/login/google/authorized",
)
app.register_blueprint(blueprint, url_prefix="/login")
blueprint.backend = SQLAlchemyBackend(OAuth, db.session, user=current_user, user_required=False)



app.register_blueprint(blueprint, url_prefix="/login")

db.create_all()

# Load user
@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))

# Routes
@app.route('/')
def home():
    return render_template('home.html')

@app.route('/login/google')
def google_login():
    if not google.authorized:
        return redirect(url_for("google.login"))
    resp = google.get("/oauth2/v1/userinfo")
    if resp.ok:
        user_data = resp.json()
        email = user_data["email"]
        name = user_data["name"]
        user = User.query.filter_by(email=email).first()
        if not user:
            user = User(name=name, email=email)
            db.session.add(user)
            db.session.commit()
        login_user(user)
        flash("You are now logged in!", "success")
        return redirect(url_for('home'))
    else:
        flash("Failed to log in with Google.", "error")
        return redirect(url_for('home'))

@app.route('/logout')
@login_required
def logout():
    logout_user()
    flash("You have been logged out.", "success")
    return redirect(url_for('home'))

@app.route('/ask', methods=['GET', 'POST'])
@login_required
def ask():
    if request.method == 'POST':
        data = request.form['data']
        question = Question(data=data, user_id=current_user.id)
        db.session.add(question)
        db.session.commit()
        flash("Your question has been submitted!", "success")
        return redirect(url_for('home'))
    return render_template('ask.html')

@app.route('/questions')
@login_required
def questions():
    questions = Question.query.filter_by(user_id=current_user.id).all()
    return render_template('questions.html', questions=questions)

@app.route('/questions/edit/<int:question_id>', methods=['GET', 'POST'])
@login_required
def edit_question(question_id):
    question = Question.query.get_or_404(question_id)
    if request.method == 'POST':
        question.data = request.form['data']
        question.response = request.form['response']
        db.session.commit()
        flash("The question has been updated!", "success")
        return redirect(url_for('questions'))
    return render_template('edit_question.html', question=question)

if __name__ == '__main__':
    app.run(debug=True)

