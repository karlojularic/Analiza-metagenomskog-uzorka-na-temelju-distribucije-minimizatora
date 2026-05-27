import pandas as pd
import os
from sklearn.metrics import classification_report, confusion_matrix

def main():
    data = pd.read_csv('output/classification_data.csv')
    data = data.dropna(subset=['Most Similar Reference', 'Actual Reference'])
    print(data.head())
    columns = ['Actual Reference', 'Most Similar Reference']

    mapping = {
        "NC_004722.1": "bacillus_cereus",
        "NC_004721.2": "bacillus_cereus",
        "NC_000913.3": "escherichia_coli",
        "NC_000907.1": "haemophilus_influenzae",
        "NC_000915.1": "helicobacter_pylori",
        "NC_008530.1": "lactobacillus_gasseri"
    }

    y_true = data[columns[0]]
    y_pred = [mapping.get(x.strip()) for x in data[columns[1]]]
    print(classification_report(y_true, y_pred))
    print(confusion_matrix(y_true, y_pred))

    labels = list(dict.fromkeys(mapping.values()))

    cm = confusion_matrix(y_true, y_pred, labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    if os.path.exists('output/confusion_matrix.csv'):
        print("Removing existing confusion_matrix.csv file...")
        os.remove('output/confusion_matrix.csv')
    print("Saving confusion matrix to output/confusion_matrix.csv...")
    cm_df.to_csv('output/confusion_matrix.csv')

if __name__ == "__main__":
    main()