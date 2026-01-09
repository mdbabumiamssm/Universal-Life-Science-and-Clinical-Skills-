import pandas as pd
import networkx as nx

class CellCommunicator:
    """
    Framework for inferring Cell-Cell Communication (Ligand-Receptor analysis).
    Simulates the logic of CellChat/CellPhoneDB.
    """
    
    def __init__(self, lr_database: pd.DataFrame):
        """
        :param lr_database: DataFrame with columns ['Ligand', 'Receptor']
        """
        self.lr_db = lr_database
        
    def infer_interactions(self, 
                          expr_data: pd.DataFrame, 
                          cell_types: pd.Series, 
                          threshold: float = 0.5) -> pd.DataFrame:
        """
        Predicts communication between cell type A and cell type B.
        Logic:
        - Ligand must be expressed in Sender (A)
        - Receptor must be expressed in Receiver (B)
        """
        cell_groups = expr_data.groupby(cell_types).mean()
        interactions = []
        
        print("Inferring Cell-Cell Interactions...")
        
        for sender in cell_groups.index:
            for receiver in cell_groups.index:
                if sender == receiver: continue # Simplified: ignore autocrine for now
                
                for _, row in self.lr_db.iterrows():
                    ligand = row['Ligand']
                    receptor = row['Receptor']
                    
                    if ligand not in cell_groups.columns or receptor not in cell_groups.columns:
                        continue
                        
                    l_expr = cell_groups.loc[sender, ligand]
                    r_expr = cell_groups.loc[receiver, receptor]
                    
                    # Mass Action Law proxy: L * R
                    score = l_expr * r_expr
                    
                    if score > threshold:
                        interactions.append({
                            'Sender': sender,
                            'Receiver': receiver,
                            'Ligand': ligand,
                            'Receptor': receptor,
                            'Score': score
                        })
                        
        return pd.DataFrame(interactions)

    def build_communication_network(self, interactions: pd.DataFrame) -> nx.DiGraph:
        """
        Converts tabular interactions into a directed graph for visualization.
        """
        G = nx.DiGraph()
        for _, row in interactions.iterrows():
            G.add_edge(row['Sender'], row['Receiver'], 
                       interaction=f"{row['Ligand']}-{row['Receptor']}",
                       weight=row['Score'])
        return G

# --- Example Usage ---
if __name__ == "__main__":
    # Mock L-R Database
    lr_db = pd.DataFrame({
        'Ligand': ['L1', 'L2'],
        'Receptor': ['R1', 'R2']
    })
    
    # Mock Expression (Clusters x Genes)
    # Cluster A expresses L1, Cluster B expresses R1
    expr_data = pd.DataFrame({
        'L1': [0.9, 0.1, 0.1],
        'R1': [0.1, 0.9, 0.1],
        'L2': [0.1, 0.1, 0.8],
        'R2': [0.1, 0.1, 0.8]
    }, index=['Cell1', 'Cell2', 'Cell3'])
    
    labels = pd.Series(['ClusterA', 'ClusterB', 'ClusterC'], index=['Cell1', 'Cell2', 'Cell3'])
    
    comm = CellCommunicator(lr_db)
    res = comm.infer_interactions(expr_data, labels, threshold=0.1)
    
    print("Significant Interactions:")
    print(res)
