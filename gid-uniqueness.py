#! /usr/bin/env python3

import csv
from collections import defaultdict

# Initialize structures to count occurrences and store units
gid_count = defaultdict(int)
gid_units = defaultdict(list)
MEMBERS_ONLY = True

# Read the CSV file
with open('/mnt/d/GID-Inventory-Export.csv', 'r', encoding='utf-8-sig') as file:  # Handle potential BOM
    reader = csv.DictReader(file)
    for row in reader:
        gid_number = row['gidNumber']
        unit = row['Unit']
        group = row['group']
        members = row['members']
        # Increment count for gidNumber
        gid_count[gid_number] += 1
        # Add the unit with group name and member count to the list for gidNumber if MEMBERS_ONLY is True and members exist
        if not MEMBERS_ONLY or members.strip():
            member_count = len(members.split(',')) if members.strip() else 0
            gid_units[gid_number].append((unit, group, member_count))

# Output the results
print("GIDNumber,Occurrences,Units,Change")
for gid, count in gid_count.items():
    if len(gid_units[gid]) > 1:  # Only print lines with more than one unit affected
        # Sort units by member count (descending) and alphabetically as a tie-breaker
        sorted_units = sorted(gid_units[gid], key=lambda x: (-x[2], x[0]))
        units_with_details = "/".join([f"{unit}({group};{members})" for unit, group, members in sorted_units])
        
        # Determine the 'Change' column value
        max_members = sorted_units[0][2]  # Largest member count
        if sum(1 for _, _, members in sorted_units if members == max_members) > 1:  # Check for ties
            change = "/".join([unit for unit, _, _ in sorted_units])  # All units if tie
        else:
            change = "/".join([unit for unit, _, _ in sorted_units if unit != sorted_units[0][0]])  # Exclude the max unit
        
        print(f"{gid},{count},{units_with_details},{change}")

