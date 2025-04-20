from utils.classes import *
import utils.check_files as cf
import os
import xgboost
import pickle
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

def classify_points(protein, feature_vector): 
    print("Started classifying points.")
    script_dir = os.path.dirname(__file__)
    model_path = os.path.join(script_dir, "xgboost_binding_site_model_20250414_115052.pkl")
    
    with open(model_path, "rb") as file:
        model = pickle.load(file)
    file.close()

    predictions = model.predict_proba(feature_vector)
    clustering_points = []

    for point, prediction in zip(protein.get_points(), predictions):
        prob = prediction[1]
        point.set_probability(prob)
        if prob >= 0.75:
            clustering_points.append(point)

    print("Finished classifying points.")

    return clustering_points

def clustering(clustering_points, distance = 3):
    clustering_coordinates = np.array([point.get_coord() for point in clustering_points])
    clustering_distances = pdist(clustering_coordinates)
    clustering =  linkage(clustering_distances, method = "single")
    cluster_labels = fcluster(clustering, t = distance, criterion = "distance")

    clusters = {}
    for point, label in zip(clustering_points, cluster_labels):
        if label not in clusters:
            clusters[label] = Cluster()

        score = point.get_probability() ** 2

        clusters[label].append_point(point)
        clusters[label].add_score(score)

    return clusters

def filtering(final_clusters, size = 3, atoms = 160):
    final_clusters = {label: cluster for label, cluster in final_clusters.items() if len(cluster.get_points()) >= size}
    labels_to_remove = []
    
    for label, cluster in final_clusters.items():
        unique_atoms = set()
        for point in cluster.get_points():
            for atom in point.get_atoms():
                unique_atoms.add(atom)
        if len(unique_atoms) > atoms:
            labels_to_remove.append(label)
        else:
            cluster.add_atoms(unique_atoms)

    for label in labels_to_remove:
        final_clusters.pop(label)

    sorted_clusters = sorted(final_clusters.items(), key=lambda item: item[1].get_score(), reverse=True)

    print("Finished clustering.")

    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "results.log")

    with open(output_file, "w") as file:
        file.write(f"There are a total of {len(sorted_clusters)} clusters with distance {size} and maximum number of atoms {atoms}\n")
        for label, cluster in sorted_clusters:
            file.write(f"Cluster {label}: Score = {cluster.get_score():.4f}\n")
    file.close()
    
    return sorted_clusters



def cluster_points(protein, feature_vector):
    clustering_points = classify_points(protein, feature_vector)
    print("Started clustering.")
    clusters = clustering(clustering_points)
    filtered_clusters = filtering(clusters)
    return filtered_clusters

