import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

# Function to generate synthetic ECG waveform
def generate_ecg_waveform(length, pattern="normal"):
    t = np.linspace(0, 1, length)
    
    if pattern == "normal":
        # Normal sinus rhythm: simple sine wave approximation
        ecg = np.sin(2 * np.pi * 5 * t) + 0.1 * np.random.normal(0, 1, length)
    elif pattern == "afib":
        # AFib: irregular waveform with noise
        ecg = np.sin(2 * np.pi * 7 * t) + 0.5 * np.random.normal(0, 1, length)
        ecg[np.random.randint(0, length, size=length//20)] += np.random.normal(0, 1, length//20)
    else:
        raise ValueError("Invalid pattern type. Choose 'normal' or 'afib'.")
    
    return t, ecg

# Function to generate random activity labels
def generate_activity_labels(num_segments):
    activities = ["stress", "caffeine", "rest", "exercise"]
    return [random.choice(activities) for _ in range(num_segments)]

# Generate ECG data
num_segments = 100  # Number of segments
segment_length = 1000  # Number of samples per segment

data = []
for i in range(num_segments):
    pattern = "normal" if i % 2 == 0 else "afib"  # Alternate between normal and afib
    t, ecg = generate_ecg_waveform(segment_length, pattern=pattern)
    activity = random.choice(["stress", "caffeine", "rest", "exercise"])
    for j in range(segment_length):
        data.append({
            "time": t[j] + i,  # Append segment index to time
            "ecg_value": ecg[j],
            "pattern": pattern,
            "activity": activity
        })

# Convert to DataFrame
df = pd.DataFrame(data)

# Save data to CSV
df.to_csv("synthetic_ecg_data.csv", index=False)
print("Synthetic ECG data saved to 'synthetic_ecg_data.csv'")