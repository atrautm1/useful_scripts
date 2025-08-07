"""
This code prioritizes typability of passwords with some security.
Useful for people who don't like to use password managers.
"""

import secrets
import string

word_list = []
number_of_words = 2
position = 1

# Step 1: Generate the base password using random words
with open("/usr/share/dict/words") as f:
    words = [word.strip() for word in f]
    word_list = [secrets.choice(words) for i in range(0, number_of_words)]

# Make the last word title case
word_list[-1] = word_list[-1].title()

# Step 2: Generate a random number and a special character
random_number = str(secrets.choice(range(10)))  # Random digit as a string
special_characters = string.punctuation  # List of special characters
random_special_char = secrets.choice(special_characters)
word_list.insert(
    position, f"{random_number}{random_special_char}"
)  # insert into list at position

# Output final password
final_password = "".join(word_list)
print(final_password)
