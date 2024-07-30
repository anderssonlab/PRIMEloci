import os
import pickle
import argparse
import numpy as np
import json
from sklearn.utils import shuffle
from lightgbm import LGBMClassifier
from python.extraction import extract_filenames
from python.profile_for_training import profile_for_training

# Set a random seed for reproducibility
seed = 0
np.random.seed(seed)

# Load configuration from JSON file
with open('model_config.json', 'r') as f:
    config = json.load(f)

# Argument parsing
parser = argparse.ArgumentParser(description='Train a LightGBM model.')
parser.add_argument('--model_name', type=str, required=True,
                    help='Name of the model to save.')
parser.add_argument('--model_dir', default='./', type=str, required=True,
                    help='Directory to save the model.')
parser.add_argument('--profile_main_dir', type=str, required=True,
                    help='Main directory for profiles.')
parser.add_argument('--keep_word', type=str, default=None,
                    help=('Comma-separated list of words to keep '
                          'in the file names.'))
parser.add_argument('--drop_word', type=str, default=None,
                    help='Comma-separated list of words to drop '
                         'from the file names.')
args = parser.parse_args()


# Command-line arguments
model_name = args.model_name
model_dir = args.model_dir
profile_main_dir = args.profile_main_dir

# Convert comma-separated strings to lists of strings
keep_word = args.keep_word.split(',') if args.keep_word else None
drop_word = args.drop_word.split(',') if args.drop_word else None

# Directory paths
pos_profile_dir = os.path.join(profile_main_dir, "profiles_subtnorm", "pos")
neg_profile_dir = os.path.join(profile_main_dir, "profiles_subtnorm", "neg")

# Extract file names
pos_profile_file_ls = extract_filenames(pos_profile_dir,
                                        keep_word=keep_word,
                                        drop_word=drop_word)
neg_profile_file_ls = extract_filenames(neg_profile_dir,
                                        keep_word=keep_word,
                                        drop_word=drop_word)

# Extract profiles and labels
pos_profile_np, pos_label_np = profile_for_training(pos_profile_dir,
                                                    pos_profile_file_ls,
                                                    1)
neg_profile_np, neg_label_np = profile_for_training(neg_profile_dir,
                                                    neg_profile_file_ls,
                                                    0)

# Combine positive and negative profiles and labels
x = np.concatenate((pos_profile_np, neg_profile_np))
y = np.concatenate((pos_label_np, neg_label_np))

# Shuffle the data
X_train, y_train = shuffle(x, y)

# Print the shape of the training data
print('Training data shape:')
print(X_train.shape)
print(y_train.shape)

# Load hyperparameters from the config dictionary
n_estimators = config.get('N_ESTIMATORS', 100)
learning_rate = config.get('LEARNING_RATE', 0.1)
num_leaves = config.get('NUM_LEAVES', 31)
max_depth = config.get('MAX_DEPTH', -1)
reg_alpha = config.get('REG_ALPHA', 0.0)
reg_lambda = config.get('REG_LAMBDA', 0.0)
min_split_gain = config.get('MIN_SPLIT_GAIN', 0.0)
boosting_type = config.get('BOOSTING_TYPE', 'gbdt')


# Fit model on train set
model = LGBMClassifier(n_estimators=n_estimators,
                       learning_rate=learning_rate,
                       num_leaves=num_leaves,
                       max_depth=max_depth,
                       reg_alpha=reg_alpha,
                       reg_lambda=reg_lambda,
                       min_split_gain=min_split_gain,
                       boosting_type=boosting_type)
model.fit(X_train, y_train)
print('Training set score: {:.4f}'.format(model.score(X_train, y_train)))

# Save the model to disk
filename = os.path.join(model_dir, model_name + '.sav')
pickle.dump(model, open(filename, 'wb'))
