#Experimental script, trying to combine data and exclude some 

#All of these are needed
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sb

#Importing data, 2b data now
path = r"C:\Users\chris\OneDrive\Dokumenter\Universitet\Master Project\Python\Data sets\plate_2b.csv"
#We find the path by shift clicking on the file and selecting "copy as pathway"
data = pd.read_csv(path)

#Defining names of each column
count_parental = data["Count_gfp_objects"]
count_resistant = data["Count_rfp_objects"]
rowletter_and_number = data["FileName_gfp"]

#Remember time is inferred here since we just have a number of images. 4h pr image
time = data["ImageNumber"]

#Using a loop, i can create a dictionary which stores each of the different letter 
#Codes as dfs[B2] and so on which can be called later.
#Can create the same for other letters by simply copying the below code and changing
#the letter from B to C and so on.
fig, axs = plt.subplots(6, 9, figsize=(50, 20))
eps = 1e-10
dfs = {}
model = {}
modelp = {}
prefixes = ["B", "C", "D", "E", "F", "G"]
for prefix in prefixes:
    for i in range(2, 11):  # Can adjust the range according to my needs
        #Here we designate a subset, meaning we say that the data that starts with the prefix
        #will be selected
        prefix_i = f'{prefix}{i}'
        subset = data[rowletter_and_number.str.startswith(prefix_i)].copy()
        
        #Here we specify that the gfp and rfp objects should be taken the
        #logarithm of and we add eps to avoid problems with log0
        subset["Count_gfp_objects"] = np.log(subset['Count_gfp_objects'] + eps)
        subset["Count_rfp_objects"] = np.log(subset['Count_rfp_objects'] + eps)
        subset["ImageNumberWithinGroup"] = range(len(subset))  # add new column that counts up from 0
        subset["ImageNumberWithinGroup"] = (subset["ImageNumberWithinGroup"] * 4)  # adjust for the 4h per image
        dfs[prefix_i] = subset
        #A for loop can also be used to create the linear regression models like the ones
        #in my previous attempt, but this creates them as a library which can be 
        #called later without having to brute force them like previously.
        #For Resistant
        modelnumber = f'{prefix}{i}'
        model[modelnumber] = LinearRegression()
        model[modelnumber].fit(dfs[prefix_i]["ImageNumberWithinGroup"].values.reshape(-1, 1), dfs[prefix_i]["Count_rfp_objects"])
        #For parenal
        modelp[modelnumber] = LinearRegression()
        modelp[modelnumber].fit(dfs[prefix_i]["ImageNumberWithinGroup"].values.reshape(-1, 1), dfs[prefix_i]["Count_gfp_objects"])

#Defining a function to annotate each graph with the linear regression line 
def plot_annotate(model, r_sq, data, x, y, label, ax, ypos_factor, xpos_factor):
    # Annotate the equation
    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()
    y_pos = y_max - ypos_factor * (y_max - y_min)  # adjust position based on ypos_factor
    x_pos = x_max - xpos_factor * (x_max - x_min)
    ax.text(x_pos, y_pos, f"{label} = {model.intercept_:.2f} + {model.coef_[0]:.2f}*Time, R² = {r_sq:.2f} ", fontsize = 8)

#Trying to create a for loop to avoid having to brute force the graphs
for j, prefix in enumerate(prefixes):
    for i in range (2, 11):
        #Making a "modelnumber" variable to store the prefix and i in to compact the code
        modelnumber = f'{prefix}{i}'
        #Using "get" to more safely retrieve the data from the library. If there is none, it will return "none"
        df = dfs.get(modelnumber)
        #Adding the below if statement to make sure the code continiues to work even if there is no data for a particular point
        if df is None or df.empty:
            continue
        row = j  # row index is just j, the index of the prefix in the list
        col = i - 2  # subtract 2 so that i=2 corresponds to col=0
        # Setting the title of the subplot to the model number
        axs[row, col].set_title(modelnumber)
        
        sb.scatterplot(data=df, x="ImageNumberWithinGroup", y="Count_gfp_objects" , color="blue", ax=axs[row, col])
        sb.scatterplot(data=df, x="ImageNumberWithinGroup", y='Count_rfp_objects', color="orange", ax=axs[row, col])
        sb.lineplot(data=df, x="ImageNumberWithinGroup", y= "Count_gfp_objects", ax=axs[row, col])
        sb.lineplot(data=df, x="ImageNumberWithinGroup", y='Count_rfp_objects', ax=axs[row, col])
        # plot the fitted line using the coordinates in the axs field first
        axs[row, col].plot(df["ImageNumberWithinGroup"], model[modelnumber].predict(df["ImageNumberWithinGroup"].values.reshape(-1, 1)), color='red')
        axs[row, col].plot(df["ImageNumberWithinGroup"], modelp[modelnumber].predict(df["ImageNumberWithinGroup"].values.reshape(-1, 1)), color='purple')
        #Plotting labels for axis and adding legend which displays the predecided labels for each graph
        axs[row, col].set_xlabel("Time(H)")
        axs[row, col].set_ylabel("Log Counts")
        #axs[row, col].legend()
        
        #Calculating R^2 value
        r_sq = model[modelnumber].score(df["ImageNumberWithinGroup"].values.reshape(-1, 1), df['Count_rfp_objects'])
        r_sqp = modelp[modelnumber].score(df["ImageNumberWithinGroup"].values.reshape(-1, 1), df["Count_gfp_objects"])
        #Adding the annotations for the linear regression lines
        plot_annotate(model[modelnumber], r_sq, df, "ImageNumberWithinGroup", 'Count_rfp_objects', "CRES", axs[row, col], 0.1, 0.96)
        plot_annotate(modelp[modelnumber], r_sqp, df, "ImageNumberWithinGroup", "Count_gfp_objects", "CPAR", axs[row, col], 0.2, 0.96)


plt.tight_layout()
#plt.suptitle("Counts at different concentrations of Gefitnib")
plt.show()
#Making a dataframe that says which columns we are using
coef_data = pd.DataFrame(columns=["Position on 96 well plate", "Slope of Resistant Culture"])

#Trying to create a table for the slopes of the resistant cells to see how they look
for j, prefix in enumerate(prefixes):
    for i in range (2, 11):
        modelnumber = f'{prefix}{i}'
        #Using "get" to more safely retrieve the data from the library. If there is none, it will return "none"
        df = dfs.get(modelnumber)
        if df is None or df.empty:
            continue
        coef_data.loc[len(coef_data)] = [modelnumber, model[modelnumber].coef_[0]]
    
# Create a figure and axes
fig, ax = plt.subplots()

# Hide axes
ax.axis('off')

# Create a table and add it to the axes
table = plt.table(cellText=coef_data.values, 
                  colLabels=coef_data.columns, 
                  cellLoc = 'center', 
                  loc='center')
#Adding a title to the table
#fig.suptitle ("Slopes for ressitant cells")

# Scale the table
table.scale(3, 3)
#ACTIVATE THE LINE BELOW TO CREATE A EXCEL FILE WITH THE DATA
#coef_data.to_excel("Resistant_slopes_2b_plate.xlsx")

plt.show()

# Same but for the parental data
coefp_data = pd.DataFrame(columns=["Position on 96 well plate", "Slope of Parental Culture"])

for j, prefix in enumerate(prefixes):
    for i in range (2, 11):
        #Calling the dfs[modelnumber]  df to make it easier to write multiple times
        modelnumber = f'{prefix}{i}'
        #Using "get" to more safely retrieve the data from the library. If there is none, it will return "none"
        df = dfs.get(modelnumber)
        if df is None or df.empty:
            continue
        coefp_data.loc[len(coefp_data)] = [modelnumber, modelp[modelnumber].coef_[0]]
    
# Create a figure and axes
fig, ax = plt.subplots()

# Hide axes
ax.axis('off')

# Create a table and add it to the axes
table = plt.table(cellText=coefp_data.values, 
                  colLabels=coefp_data.columns, 
                  cellLoc = 'center', 
                  loc='center')
#Adding a title to the table
#fig.suptitle ("Slopes for ressitant cells")

# Scale the table
table.scale(3, 3)

#ACTIVATE THE LINE BELOW TO CREATE A EXCEL FILE WITH THE DATA
#coefp_data.to_excel("Parental_slopes_2b_plate.xlsx")

plt.show()
   

