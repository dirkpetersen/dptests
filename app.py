from flask import Flask, render_template, jsonify
import random

app = Flask(__name__)

# Emoji pairs for the memory game
EMOJI_PAIRS = ['🎮', '🎯', '🎪', '🎨', '🎭', '🎪', '🎸', '🎺', 
                '🌟', '🌈', '🌺', '🌸', '🦋', '🦄', '🐉', '🐳']

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/new-game')
def new_game():
    # Create pairs and shuffle them
    cards = EMOJI_PAIRS * 2
    random.shuffle(cards)
    return jsonify({'cards': cards})

if __name__ == '__main__':
    app.run(debug=True)
