class MemoryGame {
    constructor() {
        this.cards = [];
        this.flippedCards = [];
        this.matchedPairs = 0;
        this.moves = 0;
        this.isProcessing = false;
        this.startTime = null;
        this.timerInterval = null;
        
        this.gameBoard = document.getElementById('game-board');
        this.movesElement = document.getElementById('moves');
        this.timerElement = document.getElementById('timer');
        this.matchesElement = document.getElementById('matches');
        this.newGameBtn = document.getElementById('new-game-btn');
        this.playAgainBtn = document.getElementById('play-again-btn');
        this.winModal = document.getElementById('win-modal');
        
        this.init();
    }
    
    init() {
        this.newGameBtn.addEventListener('click', () => this.startNewGame());
        this.playAgainBtn.addEventListener('click', () => {
            this.winModal.classList.remove('show');
            this.startNewGame();
        });
        
        this.startNewGame();
    }
    
    async startNewGame() {
        // Reset game state
        this.flippedCards = [];
        this.matchedPairs = 0;
        this.moves = 0;
        this.isProcessing = false;
        this.updateStats();
        
        // Clear timer
        if (this.timerInterval) {
            clearInterval(this.timerInterval);
        }
        
        // Fetch new game data
        try {
            const response = await fetch('/api/new-game');
            const data = await response.json();
            this.cards = data.cards;
            this.renderBoard();
            this.startTimer();
        } catch (error) {
            console.error('Error starting new game:', error);
        }
    }
    
    renderBoard() {
        this.gameBoard.innerHTML = '';
        
        this.cards.forEach((emoji, index) => {
            const card = document.createElement('div');
            card.className = 'card';
            card.dataset.index = index;
            card.dataset.emoji = emoji;
            
            card.innerHTML = `
                <div class="card-front">${emoji}</div>
                <div class="card-back"></div>
            `;
            
            card.addEventListener('click', () => this.flipCard(card));
            this.gameBoard.appendChild(card);
        });
    }
    
    flipCard(card) {
        if (this.isProcessing || card.classList.contains('flipped') || card.classList.contains('matched')) {
            return;
        }
        
        card.classList.add('flipped');
        this.flippedCards.push(card);
        
        if (this.flippedCards.length === 2) {
            this.isProcessing = true;
            this.moves++;
            this.updateStats();
            
            setTimeout(() => this.checkMatch(), 800);
        }
    }
    
    checkMatch() {
        const [card1, card2] = this.flippedCards;
        const emoji1 = card1.dataset.emoji;
        const emoji2 = card2.dataset.emoji;
        
        if (emoji1 === emoji2) {
            card1.classList.add('matched');
            card2.classList.add('matched');
            this.matchedPairs++;
            this.updateStats();
            
            if (this.matchedPairs === 16) {
                setTimeout(() => this.showWinModal(), 500);
            }
        } else {
            card1.classList.remove('flipped');
            card2.classList.remove('flipped');
        }
        
        this.flippedCards = [];
        this.isProcessing = false;
    }
    
    startTimer() {
        this.startTime = Date.now();
        this.timerInterval = setInterval(() => {
            const elapsed = Math.floor((Date.now() - this.startTime) / 1000);
            const minutes = Math.floor(elapsed / 60).toString().padStart(2, '0');
            const seconds = (elapsed % 60).toString().padStart(2, '0');
            this.timerElement.textContent = `${minutes}:${seconds}`;
        }, 1000);
    }
    
    updateStats() {
        this.movesElement.textContent = this.moves;
        this.matchesElement.textContent = `${this.matchedPairs}/16`;
    }
    
    showWinModal() {
        clearInterval(this.timerInterval);
        
        document.getElementById('final-moves').textContent = this.moves;
        document.getElementById('final-time').textContent = this.timerElement.textContent;
        
        this.winModal.classList.add('show');
    }
}

// Start the game when the page loads
document.addEventListener('DOMContentLoaded', () => {
    new MemoryGame();
});
