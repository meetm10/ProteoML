from sklearn import svm
from sklearn.metrics import confusion_matrix,roc_curve
import datagen
import numpy as np
from sklearn.model_selection import train_test_split
import math
import pydot
from sklearn.tree import export_graphviz
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt

X_normal, X_abnormal, Y_normal, Y_abnormal = datagen.fetch_train_test_data()

data = X_abnormal[:2300]
data = np.concatenate((data,X_normal[:2300]))

labels = Y_abnormal[:2300]
labels = np.concatenate((labels,Y_normal[:2300]))

X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.3,random_state=78)

clf = RandomForestClassifier(n_estimators=4000, max_depth=11,random_state=0,criterion='gini',class_weight={0: 1})

clf.fit(X_train, y_train)
pred = clf.predict(X_test)
print('error:',sum(abs(pred-y_test))/len(y_test))

cm = confusion_matrix(y_test,pred)
print(cm)

true_neg = cm[0][0]
false_pos = cm[0][1]
false_neg = cm[1][0]
true_pos = cm[1][1]

sensitivity = true_pos/(true_pos+false_neg)
print(sensitivity)
specificity = true_neg/(true_neg+false_pos)
print(specificity)
likelihood_ratio = sensitivity/(1-specificity)
print(likelihood_ratio)
balance_accuracy = (sensitivity+specificity)/2
print(balance_accuracy)
mcc = ((true_neg*true_pos)-(false_neg*false_pos))/math.sqrt((true_pos+false_pos)*(true_pos+false_neg)*(true_neg+false_pos)*(true_neg+false_neg))
print(mcc)

y_pred_rf = clf.predict_proba(X_test)[:, 1]
fpr_rf, tpr_rf, _ = roc_curve(y_test, y_pred_rf)

plt.title('ROC Curve')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.plot(fpr_rf, tpr_rf,'b')
plt.show()

estimator = clf.estimators_[5]

export_graphviz(estimator, out_file='tree.dot', 
                feature_names = ['Aliphatic','Aromatic','Both_charged','Hydrophobic','Neg_charge_density','Neg_charged','Neg_skew','Pos_charged','Pos_skew','Special','Uncharged'],
                class_names = ['Normal','IDP/IDR'],
                rounded = True, proportion = False, 
                precision = 2, filled = True)

(graph,) = pydot.graph_from_dot_file('tree.dot')
graph.write_png('tree.png')