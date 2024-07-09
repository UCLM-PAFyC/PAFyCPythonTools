import numpy as np

# matrix representation of the data - has shape (n, d), d are the features, n are the objects
data = np.array([
    [   1,   2,  -1,   4,  10],
    [   3,  -3,  -3,  12, -15],
    [   2,   1,  -2,   4,   5],
    [   5,   1,  -5,  10,   5],
    [   2,   3,  -3,   5,  12],
    [   4,   0,  -3,  16,   2],
])
# use `ddof = 1` if using sample data (default assumption) and use `ddof = 0` if using population data
standardized_data = (data - data.mean(axis = 0)) / data.std(axis = 0)
covariance_matrix = np.cov(standardized_data, ddof = 1, rowvar = False)
eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
# np.argsort can only provide lowest to highest; use [::-1] to reverse the list
order_of_importance = np.argsort(eigenvalues)[::-1]
# utilize the sort order to sort eigenvalues and eigenvectors
sorted_eigenvalues = eigenvalues[order_of_importance]
sorted_eigenvectors = eigenvectors[:,order_of_importance] # sort the columns
# use sorted_eigenvalues to ensure the explained variances correspond to the eigenvectors
explained_variance = sorted_eigenvalues / np.sum(sorted_eigenvalues)
k = 2 # select the number of principal components
reduced_data = np.matmul(standardized_data, sorted_eigenvectors[:,:k]) # transform the original data
total_explained_variance = sum(explained_variance[:k])
yo = 1

