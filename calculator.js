let currentInput = '0';
let history = '';
let result = null;
let lastOperation = '';
let decimalAdded = false;

// DOM elements
const currentDisplay = document.getElementById('current');
const historyDisplay = document.getElementById('history');

// Update the display
function updateDisplay() {
    currentDisplay.textContent = currentInput;
    historyDisplay.textContent = history;
}

// Append a number to the current input
function appendNumber(number) {
    if (currentInput === '0' || currentInput === 'Error') {
        currentInput = number;
    } else {
        currentInput += number;
    }
    updateDisplay();
}

// Add a decimal point
function appendDecimal() {
    if (!decimalAdded) {
        currentInput += '.';
        decimalAdded = true;
        updateDisplay();
    }
}

// Handle operators
function appendOperator(operator) {
    // Special case for square root
    if (operator === 'Math.sqrt') {
        history = `√(${currentInput})`;
        try {
            currentInput = String(Math.sqrt(parseFloat(currentInput)));
            if (currentInput === 'NaN') {
                currentInput = 'Error';
            }
        } catch (e) {
            currentInput = 'Error';
        }
        decimalAdded = currentInput.includes('.');
        updateDisplay();
        return;
    }

    // For other operators
    if (currentInput !== 'Error') {
        // If we have a pending operation, calculate it first
        if (result !== null) {
            calculate();
        } else {
            result = parseFloat(currentInput);
        }
        
        history = currentInput + ' ' + operator;
        lastOperation = operator;
        currentInput = '0';
        decimalAdded = false;
        updateDisplay();
    }
}

// Calculate the result
function calculate() {
    if (currentInput === 'Error' || lastOperation === '') {
        return;
    }

    let currentValue = parseFloat(currentInput);
    history += ' ' + currentInput;

    try {
        switch (lastOperation) {
            case '+':
                result += currentValue;
                break;
            case '-':
                result -= currentValue;
                break;
            case '*':
                result *= currentValue;
                break;
            case '/':
                if (currentValue === 0) {
                    throw new Error('Division by zero');
                }
                result /= currentValue;
                break;
            case '%':
                result %= currentValue;
                break;
        }

        currentInput = String(result);
        if (currentInput === 'NaN' || currentInput === 'Infinity') {
            currentInput = 'Error';
        }
    } catch (e) {
        currentInput = 'Error';
    }

    // Reset for next calculation
    lastOperation = '';
    decimalAdded = currentInput.includes('.');
    updateDisplay();
}

// Clear all data
function clearAll() {
    currentInput = '0';
    history = '';
    result = null;
    lastOperation = '';
    decimalAdded = false;
    updateDisplay();
}

// Delete the last character
function deleteChar() {
    if (currentInput === 'Error') {
        clearAll();
        return;
    }
    
    if (currentInput.length === 1) {
        currentInput = '0';
    } else {
        if (currentInput[currentInput.length - 1] === '.') {
            decimalAdded = false;
        }
        currentInput = currentInput.slice(0, -1);
    }
    updateDisplay();
}

// Initialize display
updateDisplay();
