def generate_EQXGD_data(n_samples=1000, n_per_sample=50):
    """
    Generates data for training an ANN model based on the EQXGD distribution.

    Parameters:
    - n_samples: Number of sample sets to generate
    - n_per_sample: Number of samples in each set

    Returns:
    - X: r_samples (the input data)
    - Y: Theta, Alpha, Beta (the output parameters)
    """
    # Randomly generate theta, alpha, beta values for each sample from uniform distributions
    theta_values = np.random.uniform(0.1, 5, n_samples)  # Example range for theta
    alpha_values = np.random.uniform(0.1, 5, n_samples)  # Example range for alpha
    beta_values = np.random.uniform(0.1, 5, n_samples)   # Example range for beta
    # Generate r_samples using the generated theta, alpha, beta values
    r_samples = np.zeros((n_samples, n_per_sample))
    for i in range(n_samples):
        r_samples[i, :] = rEQXGD(n_per_sample, theta_values[i], alpha_values[i], beta_values[i])
    # Combine the data: r_sample will be the input (X), theta, alpha, beta will be the output (Y)
    X = r_samples
    X.sort(axis=1)
    Y = np.vstack([theta_values, alpha_values, beta_values]).T
    return X, Y
X,Y=generate_EQXGD_data(n_samples=10000, n_per_sample=30)



import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Input
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

def train_ann_model(X, Y, model_save_path='ann_model.h5', epochs=2000, batch_size=64, test_size=0.2, random_state=123):
    """
    Train an Artificial Neural Network (ANN) model to predict parameters based on input data.

    Parameters:
    X (numpy.ndarray): Input features, a 2D array where each row represents a sample.
    Y (numpy.ndarray): Target output, a 2D array where each row corresponds to the parameters (theta, alpha, beta) for each sample.
    model_save_path (str): Path to save the trained model. Default is 'ann_model.h5'.
    epochs (int): Number of training epochs. Default is 50.
    batch_size (int): Number of samples per gradient update. Default is 64.
    test_size (float): Proportion of the dataset to include in the test split. Default is 0.2.
    random_state (int): Seed for the random number generator. Default is 123.

    Returns:
    model: The trained ANN model.
    history: Training history containing loss values and metrics over epochs.
    predictions: Predictions made by the model on the test set.
    actual: Actual values of theta, alpha, and beta from the test set.

    Example usage:
    model, history, predictions, actual = train_ann_model(X, Y)
    """

    # Split the data into training and test sets
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=random_state)

    # Normalize the input features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Build the ANN model
    model = Sequential()
    model.add(Input(shape=(X.shape[1],)))  # Explicit Input layer
    model.add(Dense(128, activation='relu'))  # First hidden layer with 128 neurons
    model.add(Dense(64, activation='relu'))   # Second hidden layer with 64 neurons
    model.add(Dense(32, activation='relu'))   # Third hidden layer with 32 neurons
    model.add(Dense(3, activation='relu'))     # Output layer: 3 neurons for theta, alpha, and beta

    # Compile the model with Adam optimizer and mean squared error loss
    model.compile(optimizer='adam', loss='mean_squared_error')

    # Train the model on the training data
    history = model.fit(X_train, Y_train, epochs=epochs, batch_size=batch_size, validation_split=0.2)

    # Evaluate the model on the test set
    test_loss = model.evaluate(X_test, Y_test)
    print(f"Test loss: {test_loss}")

    # Make predictions on the test set
    predictions = model.predict(X_test)
    print(f"Predicted theta, alpha, beta: \n{predictions[:5]}")
    print(f"Actual theta, alpha, beta: \n{Y_test[:5]}")

    # Save the trained model to the specified path
    model.save(model_save_path)
    print(f"Model saved to {model_save_path}")

    return model, history, predictions, Y_test
