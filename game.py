import random

def number_guessing_game():
    """
    A simple number guessing game where the player tries to guess
    a random number between 1 and 100.
    """
    print("Welcome to the Number Guessing Game!")
    print("I'm thinking of a number between 1 and 100.")
    
    # Generate a random number between 1 and 100
    secret_number = random.randint(1, 100)
    attempts = 0
    max_attempts = 10
    
    while attempts < max_attempts:
        # Get player's guess
        try:
            guess = int(input(f"Attempt {attempts+1}/{max_attempts}. Enter your guess: "))
        except ValueError:
            print("Please enter a valid number.")
            continue
        
        attempts += 1
        
        # Check the guess
        if guess < secret_number:
            print("Too low!")
        elif guess > secret_number:
            print("Too high!")
        else:
            print(f"Congratulations! You guessed the number in {attempts} attempts!")
            break
    
    if attempts == max_attempts and guess != secret_number:
        print(f"Game over! You've used all {max_attempts} attempts.")
        print(f"The secret number was {secret_number}.")

if __name__ == "__main__":
    play_again = "y"
    while play_again.lower() == "y":
        number_guessing_game()
        play_again = input("Do you want to play again? (y/n): ")
    
    print("Thanks for playing!")
