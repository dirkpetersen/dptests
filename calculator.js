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
let currentInput = '0';
let previousInput = '';
let operation = null;
let resetScreen = false;

function updateDisplay() {
    document.getElementById('display').textContent = currentInput;
}

function appendToDisplay(value) {
    if (currentInput === '0' || resetScreen) {
        currentInput = value;
        resetScreen = false;
    } else {
        // Prevent multiple decimal points
        if (value === '.' && currentInput.includes('.')) return;
        
        // Limit display length to prevent overflow
        if (currentInput.length < 12) {
            currentInput += value;
        }
    }
    updateDisplay();
}

function clearDisplay() {
    currentInput = '0';
    previousInput = '';
    operation = null;
    updateDisplay();
}

function backspace() {
    if (currentInput.length > 1) {
        currentInput = currentInput.slice(0, -1);
    } else {
        currentInput = '0';
    }
    updateDisplay();
}

function toggleSign() {
    currentInput = (parseFloat(currentInput) * -1).toString();
    updateDisplay();
}

function handleOperator(op) {
    const value = parseFloat(currentInput);
    
    if (previousInput === '') {
        previousInput = currentInput;
    } else if (operation) {
        const result = performCalculation();
        previousInput = result.toString();
        addToHistory(`${previousInput} ${op}`);
    }
    
    operation = op;
    resetScreen = true;
}

function calculate() {
    if (!operation || previousInput === '') return;
    
    const result = performCalculation();
    const calculation = `${previousInput} ${operation} ${currentInput} = ${result}`;
    
    currentInput = result.toString();
    previousInput = '';
    operation = null;
    
    addToHistory(calculation);
    updateDisplay();
}

function performCalculation() {
    const prev = parseFloat(previousInput);
    const current = parseFloat(currentInput);
    
    let result;
    switch (operation) {
        case '+':
            result = prev + current;
            break;
        case '-':
            result = prev - current;
            break;
        case '*':
            result = prev * current;
            break;
        case '/':
            result = prev / current;
            break;
        case '%':
            result = prev % current;
            break;
        default:
            return current;
    }
    
    // Handle precision issues
    return Math.round(result * 1000000) / 1000000;
}

function addToHistory(calculation) {
    const historyDiv = document.getElementById('history');
    const historyItem = document.createElement('div');
    historyItem.className = 'history-item';
    historyItem.textContent = calculation;
    historyDiv.appendChild(historyItem);
    
    // Scroll to bottom of history
    historyDiv.scrollTop = historyDiv.scrollHeight;
}

// Initialize
updateDisplay();
