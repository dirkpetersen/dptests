let selectedUsers = new Set();
let selectedGroup = null;

function searchUsers() {
    const query = document.getElementById('userSearchInput').value;
    fetch(`/api/users/search?q=${encodeURIComponent(query)}`)
        .then(response => response.json())
        .then(users => {
            const resultsDiv = document.getElementById('userSearchResults');
            resultsDiv.innerHTML = users.map(user => `
                <div class="user-card ${selectedUsers.has(user.id) ? 'selected' : ''}" 
                     onclick="toggleUser(${user.id}, '${user.name}', '${user.email}')">
                    <h5>${user.name}</h5>
                    <p>${user.email}</p>
                    <small>${user.department}</small>
                </div>
            `).join('');
        });
}

function searchGroups() {
    const query = document.getElementById('groupSearchInput').value;
    fetch(`/api/groups/search?q=${encodeURIComponent(query)}`)
        .then(response => response.json())
        .then(groups => {
            const resultsDiv = document.getElementById('groupSearchResults');
            resultsDiv.innerHTML = groups.map(group => `
                <div class="group-card ${selectedGroup?.id === group.id ? 'selected' : ''}"
                     onclick="selectGroup(${group.id}, '${group.name}')">
                    <h5>${group.name}</h5>
                    <p>${group.description}</p>
                </div>
            `).join('');
        });
}

function toggleUser(id, name, email) {
    if (selectedUsers.has(id)) {
        selectedUsers.delete(id);
    } else {
        selectedUsers.add(id);
    }
    updateSelectedUsers();
    document.querySelectorAll(`.user-card`).forEach(card => {
        if (card.onclick.toString().includes(`${id}`)) {
            card.classList.toggle('selected');
        }
    });
}

function selectGroup(id, name) {
    selectedGroup = { id, name };
    updateSelectedGroup();
    document.querySelectorAll('.group-card').forEach(card => {
        card.classList.remove('selected');
        if (card.onclick.toString().includes(`${id}`)) {
            card.classList.add('selected');
        }
    });
}

function updateSelectedUsers() {
    const container = document.getElementById('selectedUsers');
    container.innerHTML = Array.from(selectedUsers).map(id => {
        const userCard = document.querySelector(`.user-card[onclick*="${id}"]`);
        const name = userCard.querySelector('h5').textContent;
        return `
            <span class="selected-item">
                ${name}
                <span class="remove-btn" onclick="toggleUser(${id})">×</span>
            </span>
        `;
    }).join('');
}

function updateSelectedGroup() {
    const container = document.getElementById('selectedGroup');
    if (selectedGroup) {
        container.innerHTML = `
            <span class="selected-item">
                ${selectedGroup.name}
                <span class="remove-btn" onclick="selectGroup(null)">×</span>
            </span>
        `;
    } else {
        container.innerHTML = '';
    }
}

function submitChanges() {
    if (selectedUsers.size === 0 || !selectedGroup) {
        alert('Please select both users and a group before submitting');
        return;
    }

    const changes = {
        users: Array.from(selectedUsers),
        group: selectedGroup.id
    };

    fetch('/api/submit-changes', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify(changes)
    })
    .then(response => response.json())
    .then(result => {
        alert('Changes submitted successfully!');
        // Clear selections
        selectedUsers.clear();
        selectedGroup = null;
        updateSelectedUsers();
        updateSelectedGroup();
        // Reset search results
        document.getElementById('userSearchResults').innerHTML = '';
        document.getElementById('groupSearchResults').innerHTML = '';
    })
    .catch(error => {
        alert('Error submitting changes');
        console.error('Error:', error);
    });
}

// Initialize Bootstrap tabs
document.addEventListener('DOMContentLoaded', function() {
    var triggerTabList = [].slice.call(document.querySelectorAll('#mainTabs button'));
    triggerTabList.forEach(function(triggerEl) {
        new bootstrap.Tab(triggerEl);
    });
});
