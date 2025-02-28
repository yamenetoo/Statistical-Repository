import torch
import numpy as np
from  tqdm import tqdm
epsilon = 1e-5
def mixture_pdf(x, params, epsilon=epsilon):
    theta, a, b, w = params
    theta = torch.clamp(theta, min=epsilon)
    a = torch.clamp(a, min=epsilon)
    b = torch.clamp(b, min=epsilon)
    w = torch.clamp(w, min=epsilon)
    w1 = w / (1 + w)
    w2 = 1 / (1 + w)
    x_clamped = torch.clamp(x, min=epsilon)
    pdf_exp = theta * torch.exp(-theta * x_clamped)
    log_gamma_term = a * torch.log(b) - torch.lgamma(a) + (a - 1) * torch.log(x_clamped) - b * x_clamped
    pdf_gamma = torch.exp(log_gamma_term)
    combined_pdf = w1 * pdf_exp + w2 * pdf_gamma
    combined_pdf = torch.clamp(combined_pdf, min=epsilon)
    return combined_pdf
def negative_log_likelihood(params, data, epsilon=epsilon):
    pdf_values = mixture_pdf(data, params, epsilon)
    nll = -torch.sum(torch.log(pdf_values))
    return nll

def SGD(data, init_params, learning_rate=0.001, epochs=100, ver=False):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data_tensor = torch.tensor(data, dtype=torch.float32, device=device)
    params = torch.tensor(init_params, dtype=torch.float32, requires_grad=True, device=device)
    optimizer = torch.optim.Adam([params], lr=learning_rate)
    for epoch in range(epochs):
        optimizer.zero_grad()
        try:
            loss = negative_log_likelihood(params, data_tensor)
            if torch.isnan(loss) or torch.isinf(loss):
                print("Loss became NaN or Inf! Aborting.")
                return None
            loss.backward()
            torch.nn.utils.clip_grad_norm_(params, max_norm=1.0)  # Gradient clipping
            optimizer.step()
            with torch.no_grad():
                params[:] = torch.clamp(params, min=epsilon)
            if ver and epoch % 100 == 0:
                print(f"Epoch {epoch}: Loss = {loss.item():.4f}, Params = {params.tolist()}")
        except RuntimeError as e:
            print(f"Runtime error detected: {e}. Aborting.")
            return None
    return params.detach().cpu().numpy()

def generate_mixture_samples(theta, a, b, w, n_samples=1000):
    w1 = w / (1 + w)
    w2 = 1 / (1 + w)
    samples = []
    for _ in range(n_samples):
        if np.random.rand() < w1:
            samples.append(np.random.exponential(1 / theta))
        else:
            samples.append(np.random.gamma(a, 1 / b))
    return np.array(samples)


def monte_carlo_simulation(
    theta_true, a_true, b_true, w_true,
    n_simulations=100, n_samples_per_simulation=1000,
    initial_params=[0.5, 1.0, 1.5, 1.0],
    learning_rate=0.01, epochs=1000, verbose=False
):
    """
    Perform a Monte Carlo simulation to estimate parameters of a mixture model.

    Args:
        theta_true (float): True value of θ.
        a_true (float): True value of a.
        b_true (float): True value of b.
        w_true (float): True value of w.
        n_simulations (int): Number of simulations.
        n_samples_per_simulation (int): Number of samples per dataset.
        initial_params (list): Initial parameter guesses for optimization.
        learning_rate (float): Learning rate for SGD.
        epochs (int): Number of epochs for SGD.
        verbose (bool): Whether to print progress during simulations.

    Returns:
        dict: A dictionary containing the mean and standard deviation of estimated parameters.
    """
    estimated_theta = []
    estimated_a = []
    estimated_b = []
    estimated_w = []

    for i in tqdm(range(n_simulations)):
        if verbose:
            print(f"Simulation {i + 1}/{n_simulations}")
        data = generate_mixture_samples(theta_true, a_true, b_true, w_true, n_samples=n_samples_per_simulation)
        final_params = SGD(data, initial_params, learning_rate=learning_rate, epochs=epochs, ver=False)
        if final_params is not None:
            th, a_est, b_est, w_est = final_params
            estimated_theta.append(th)
            estimated_a.append(a_est)
            estimated_b.append(b_est)
            estimated_w.append(w_est)
        else:
            if verbose:
                print(f"Optimization failed for simulation {i + 1}. Skipping...")
    estimated_theta = np.array(estimated_theta)
    estimated_a = np.array(estimated_a)
    estimated_b = np.array(estimated_b)
    estimated_w = np.array(estimated_w)
    results = {
        "theta": {"mean": np.mean(estimated_theta), "std": np.std(estimated_theta)},
        "a": {"mean": np.mean(estimated_a), "std": np.std(estimated_a)},
        "b": {"mean": np.mean(estimated_b), "std": np.std(estimated_b)},
        "w": {"mean": np.mean(estimated_w), "std": np.std(estimated_w)},
    }
    return results



theta_true = 0.5
a_true = 1
b_true = 1.5
w_true = 2.5
results = monte_carlo_simulation(theta_true, a_true, b_true, w_true,n_simulations=1000, n_samples_per_simulation=1000,initial_params=[0.5, 1.0, 1.5, 1.0],learning_rate=0.01, epochs=1000, verbose=False)
print("\nMonte Carlo Simulation Results:")
print(f"True Parameters: θ={theta_true}, a={a_true}, b={b_true}, w={w_true}")
print(f"Estimated θ: Mean={results['theta']['mean']:.4f}, Std={results['theta']['std']:.4f}")
print(f"Estimated a: Mean={results['a']['mean']:.4f}, Std={results['a']['std']:.4f}")
print(f"Estimated b: Mean={results['b']['mean']:.4f}, Std={results['b']['std']:.4f}")
print(f"Estimated w: Mean={results['w']['mean']:.4f}, Std={results['w']['std']:.4f}")
