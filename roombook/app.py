from flask import Flask, request, render_template, redirect, url_for
from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///site.db'
db = SQLAlchemy(app)

class Room(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    floorplan = db.Column(db.String(120), nullable=False)
    beds = db.Column(db.Integer, nullable=False)
    bookings = db.relationship('Booking', backref='room', lazy=True)


class Booking(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    member_name = db.Column(db.String(120), nullable=False)
    guests = db.Column(db.Integer, nullable=False, default=0)
    date = db.Column(db.Date, nullable=False)
    room_id = db.Column(db.Integer, db.ForeignKey('room.id'), nullable=False)

# Initialize the Database
@app.before_first_request
def setup_db():
    db.create_all()

# Create Dummy Data
@app.before_first_request
def create_rooms():
    if not Room.query.first():
        for i in range(14):
            room = Room(floorplan=f'Floorplan {i + 1}', beds=i % 4 + 1)
            db.session.add(room)
        db.session.commit()


@app.route('/')
def index():
    rooms = Room.query.all()
    return render_template('index.html', rooms=rooms)


@app.route('/book/<int:room_id>', methods=['GET', 'POST'])
def book_room(room_id):
    room = Room.query.get_or_404(room_id)
    if request.method == 'POST':
        member_name = request.form['member_name']
        guests = int(request.form['guests'])
        date = request.form['date']
        booking = Booking(member_name=member_name, guests=guests, date=datetime.strptime(date, '%Y-%m-%d').date(), room_id=room.id)
        db.session.add(booking)
        db.session.commit()
        return redirect(url_for('index'))
    return render_template('book_room.html', room=room)


if __name__ == '__main__':
    app.run(debug=True)



