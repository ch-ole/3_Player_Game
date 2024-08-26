#Linear regression implemented, trying different way of plotting the 
#regression lines and annotating them

#Implementing logarithm transformation and linear regression for all graphs
#Lineplots with scatterplot and log transformation 
#Linear regression on top of log transformed data
#Plotting time on the x axis and 
#count of parental and resistant cells on y axis for
#ALL fractions

#All of these are needed
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sb


#We have now imported pandas as pd and can use that as our shorthand
path = r"C:\Users\chris\OneDrive\Dokumenter\Python Learning\Data sets\Raw data for learning.csv"
#We find the path by shift clicking on the file and selecting "copy as pathway"
data = pd.read_csv(path)

#Defining names of each column
fraction_parental = data["Fraction_Parental"]
replicate_id = data.iloc[:, 1]
time = data["Time"]

counts_parental = data["Counts_Parental"]

counts_resistant = data["Counts_Resistant"]

#Adding a small amount to the fractions to avoid taking the logarithm of 0
eps = 1e-10
#creating counts_parental and counts_resistant for all different fractions
#Using np.log to take the logarithm of each dataset
#Fraction 0
fraction0 = data[fraction_parental == 0]
# Multiplying by 0 here since we want the parental cell count to be 0 as there were none to begin with
fraction0["Counts_Parental"] = (np.log(fraction0["Counts_Parental"] + eps) * 0)
fraction0["Counts_Resistant"] = np.log(fraction0["Counts_Resistant"] + eps)
#Fraction 0.1
fraction10 = data[data["Fraction_Parental"] == 0.1]
fraction10["Counts_Parental"] = np.log(fraction10["Counts_Parental"]+ eps)
fraction10["Counts_Resistant"] = np.log(fraction10["Counts_Resistant"]+ eps)
#Fraction 0.3
fraction30 = data[data["Fraction_Parental"] == 0.3]
fraction30["Counts_Parental"] = np.log(fraction30["Counts_Parental"]+ eps)
fraction30["Counts_Resistant"] = np.log(fraction30["Counts_Resistant"]+ eps)
#Fraction 0.5
fraction50 = data[data["Fraction_Parental"] == 0.5]
fraction50["Counts_Parental"] = np.log(fraction50["Counts_Parental"]+ eps)
fraction50["Counts_Resistant"] = np.log(fraction50["Counts_Resistant"]+ eps)
#Fraction 0.7
fraction70 = data[data["Fraction_Parental"] == 0.7]
fraction70["Counts_Parental"] = np.log(fraction70["Counts_Parental"]+ eps)
fraction70["Counts_Resistant"] = np.log(fraction70["Counts_Resistant"]+ eps)
#Fraction 0.9
fraction90 = data[data["Fraction_Parental"] == 0.9]
fraction90["Counts_Parental"] = np.log(fraction90["Counts_Parental"]+ eps)
fraction90["Counts_Resistant"] = np.log(fraction90["Counts_Resistant"]+ eps)
#Fraction 0.95
fraction95 = data[data["Fraction_Parental"] == 0.95]
fraction95["Counts_Parental"] = np.log(fraction95["Counts_Parental"]+ eps)
fraction95["Counts_Resistant"] = np.log(fraction95["Counts_Resistant"]+ eps)
#Fraction 1
fraction1 = data[data["Fraction_Parental"] == 1]
fraction1["Counts_Parental"] = np.log(fraction1["Counts_Parental"] + eps)
#Multiplying by 0 here since there were 0 resistant cells to begin with
fraction1["Counts_Resistant"] = (np.log(fraction1["Counts_Resistant"] + eps) * 0)

#Creating an instance of the linear regression as a model
#For Parental
model10P = LinearRegression()
#Fitting the data by making a linear regression for resistant vs time in this case
model10P.fit(fraction10["Time"].values.reshape(-1, 1), fraction10["Counts_Parental"])
#For resistant
model10 = LinearRegression()
model10.fit(fraction10["Time"].values.reshape(-1, 1), fraction10["Counts_Resistant"])
# 0% Fraction
model0P = LinearRegression()
model0P.fit(fraction0["Time"].values.reshape(-1, 1), fraction0["Counts_Parental"])
model0 = LinearRegression()
model0.fit(fraction0["Time"].values.reshape(-1, 1), fraction0["Counts_Resistant"])
# 30% Fraction
model30P = LinearRegression()
model30P.fit(fraction30["Time"].values.reshape(-1, 1), fraction30["Counts_Parental"])
model30 = LinearRegression()
model30.fit(fraction30["Time"].values.reshape(-1, 1), fraction30["Counts_Resistant"])
# 50% Fraction
model50P = LinearRegression()
model50P.fit(fraction50["Time"].values.reshape(-1, 1), fraction50["Counts_Parental"])
model50 = LinearRegression()
model50.fit(fraction50["Time"].values.reshape(-1, 1), fraction50["Counts_Resistant"])
# 70% Fraction
model70P = LinearRegression()
model70P.fit(fraction70["Time"].values.reshape(-1, 1), fraction70["Counts_Parental"])
model70 = LinearRegression()
model70.fit(fraction70["Time"].values.reshape(-1, 1), fraction70["Counts_Resistant"])
# 90% Fraction
model90P = LinearRegression()
model90P.fit(fraction90["Time"].values.reshape(-1, 1), fraction90["Counts_Parental"])
model90 = LinearRegression()
model90.fit(fraction90["Time"].values.reshape(-1, 1), fraction90["Counts_Resistant"])
# 95% Fraction
model95P = LinearRegression()
model95P.fit(fraction95["Time"].values.reshape(-1, 1), fraction95["Counts_Parental"])
model95 = LinearRegression()
model95.fit(fraction95["Time"].values.reshape(-1, 1), fraction95["Counts_Resistant"])
# 100% Fraction
model1P = LinearRegression()
model1P.fit(fraction1["Time"].values.reshape(-1, 1), fraction1["Counts_Parental"])
model1 = LinearRegression()
model1.fit(fraction1["Time"].values.reshape(-1, 1), fraction1["Counts_Resistant"])

#Creating an array of 4 rows and 2 columns 
fig, axs = plt.subplots(4, 2, figsize=(10, 10))

#The code below uses 2 different plot functions. The first is "scatterplot" 
#which shows only the actual datapoints on the graph. 
#The second is "lineplot" which makes interpolated lines, as well 
#as the "shadow" which indicated where the data differs and has 
#a higher deviation from the mean

#Fraction 0
sb.scatterplot(data=fraction0, x="Time", y= "Counts_Parental", color="blue", ax=axs[0, 0])
sb.scatterplot(data=fraction0, x="Time", y="Counts_Resistant", color="orange", ax=axs[0, 0])
sb.lineplot(data=fraction0, x="Time", y= "Counts_Parental", label="CPAR 0%", ax=axs[0, 0])
sb.lineplot(data=fraction0, x="Time", y="Counts_Resistant", label="CRES 100%", ax=axs[0, 0])
# plot the fitted line using the coordinates in the axs field first
axs[0, 0].plot(fraction0["Time"], model0.predict(fraction0["Time"].values.reshape(-1, 1)), color='red')
axs[0, 0].plot(fraction0["Time"], model0P.predict(fraction0["Time"].values.reshape(-1, 1)), color='purple')
#Plotting labels for axis and adding legend which displays the predecided labels for each graph
axs[0, 0].set_xlabel("Time(H)")
axs[0, 0].set_ylabel("Log Counts")
axs[0, 0].legend()

#Fraction 0.1
sb.scatterplot(data=fraction10, x="Time", y="Counts_Parental", color="blue", ax=axs[1, 0])
sb.scatterplot(data=fraction10, x="Time", y="Counts_Resistant", color="orange", ax=axs[1, 0])
sb.lineplot(data=fraction10, x="Time", y="Counts_Parental", label="CPAR 10%", ax=axs[1, 0])
sb.lineplot(data=fraction10, x="Time", y="Counts_Resistant", label="CRES 90%", ax=axs[1, 0])
axs[1, 0].plot(fraction10["Time"], model10.predict(fraction10["Time"].values.reshape(-1, 1)), color='red')
axs[1, 0].plot(fraction10["Time"], model10P.predict(fraction10["Time"].values.reshape(-1, 1)), color='purple')
axs[1, 0].set_xlabel("Time(H)")
axs[1, 0].set_ylabel("Log Counts")
axs[1, 0].legend()

#Fraction 0.3
sb.scatterplot(data=fraction30, x="Time", y="Counts_Parental", color="blue", ax=axs[2, 0])
sb.scatterplot(data=fraction30, x="Time", y="Counts_Resistant", color="orange", ax=axs[2, 0])
sb.lineplot(fraction30, x="Time", y="Counts_Parental", label="CPAR 30%", ax=axs[2, 0])
sb.lineplot(data=fraction30, x="Time", y="Counts_Resistant", label="CRES 70%", ax=axs[2, 0])
axs[2, 0].plot(fraction30["Time"], model30.predict(fraction30["Time"].values.reshape(-1, 1)), color='red')
axs[2, 0].plot(fraction30["Time"], model30P.predict(fraction30["Time"].values.reshape(-1, 1)), color='purple')
axs[2, 0].set_xlabel("Time(H)")
axs[2, 0].set_ylabel("Log Counts")
axs[2, 0].legend()

#Fraction 0.5
sb.scatterplot(data=fraction50, x="Time", y="Counts_Parental", color="blue", ax=axs[3, 0])
sb.scatterplot(data=fraction50, x="Time", y="Counts_Resistant", color="orange", ax=axs[3, 0])
sb.lineplot(fraction50, x="Time", y="Counts_Parental", label="CPAR 50%", ax=axs[3, 0])
sb.lineplot(fraction50, x="Time", y="Counts_Resistant", label="CRES 50%", ax=axs[3, 0])
axs[3, 0].plot(fraction50["Time"], model50.predict(fraction50["Time"].values.reshape(-1, 1)), color='red')
axs[3, 0].plot(fraction50["Time"], model50P.predict(fraction50["Time"].values.reshape(-1, 1)), color='purple')
axs[3, 0].set_xlabel("Time(H)")
axs[3, 0].set_ylabel("Log Counts")
axs[3, 0].legend()

#Fraction 0.7
sb.scatterplot(data=fraction70, x="Time", y="Counts_Parental", color="blue", ax=axs[0, 1])
sb.scatterplot(data=fraction70, x="Time", y="Counts_Resistant", color="orange", ax=axs[0, 1])
sb.lineplot(fraction70, x="Time", y="Counts_Parental", label="CPAR 70%", ax=axs[0, 1])
sb.lineplot(fraction70, x="Time", y="Counts_Resistant", label="CRES 70%", ax=axs[0, 1])
axs[0, 1].plot(fraction70["Time"], model70.predict(fraction70["Time"].values.reshape(-1, 1)), color='red')
axs[0, 1].plot(fraction70["Time"], model70P.predict(fraction70["Time"].values.reshape(-1, 1)), color='purple')
axs[0, 1].set_xlabel("Time(H)")
axs[0, 1].set_ylabel("Log Counts")
axs[0, 1].legend()

#Fraction 0.9
sb.scatterplot(data=fraction90, x="Time", y="Counts_Parental", color="blue", ax=axs[1, 1])
sb.scatterplot(data=fraction90, x="Time", y="Counts_Resistant", color="orange", ax=axs[1, 1])
sb.lineplot(fraction90, x="Time", y="Counts_Parental", label="CPAR 90%", ax=axs[1, 1])
sb.lineplot(fraction90, x="Time", y="Counts_Resistant", label="CRES 10%", ax=axs[1, 1])
axs[1, 1].plot(fraction90["Time"], model90.predict(fraction90["Time"].values.reshape(-1, 1)), color='red')
axs[1, 1].plot(fraction90["Time"], model90P.predict(fraction90["Time"].values.reshape(-1, 1)), color='purple')
axs[1, 1].set_xlabel("Time(H)")
axs[1, 1].set_ylabel("Log Counts")
axs[1, 1].legend()

#Fraction 0.95
sb.scatterplot(data=fraction95, x="Time", y="Counts_Parental", color="blue", ax=axs[2, 1])
sb.scatterplot(data=fraction95, x="Time", y="Counts_Resistant", color="orange", ax=axs[2, 1])
sb.lineplot(fraction95, x="Time", y="Counts_Parental", label="CPAR 95%", ax=axs[2, 1])
sb.lineplot(fraction95, x="Time", y="Counts_Resistant", label="CRES 5%", ax=axs[2, 1])
axs[2, 1].plot(fraction95["Time"], model95.predict(fraction95["Time"].values.reshape(-1, 1)), color='red')
axs[2, 1].plot(fraction95["Time"], model95P.predict(fraction95["Time"].values.reshape(-1, 1)), color='purple')
axs[2, 1].set_xlabel("Time(H)")
axs[2, 1].set_ylabel("Log Counts")
axs[2, 1].legend()

#Fraction 1
sb.scatterplot(data=fraction1, x="Time", y="Counts_Parental", color="blue", ax=axs[3, 1])
sb.scatterplot(data=fraction1, x="Time", y= "Counts_Resistant", color="orange", ax=axs[3, 1])
sb.lineplot(fraction1, x="Time", y="Counts_Parental", label="CPAR 100%", ax=axs[3, 1])
sb.lineplot(fraction1, x="Time", y= "Counts_Resistant", label="CRES 0%", ax=axs[3, 1])
axs[3, 1].plot(fraction1["Time"], model1.predict(fraction1["Time"].values.reshape(-1, 1)), color='red')
axs[3, 1].plot(fraction1["Time"], model1P.predict(fraction1["Time"].values.reshape(-1, 1)), color='purple')
axs[3, 1].set_xlabel("Time(H)")
axs[3, 1].set_ylabel("Log Counts")
axs[3, 1].legend()

#I then want to record some values to be able to interpret the linear regression
#Starting with 0% Parental fraction
#For R^2 value:
r_sq0 = model0.score(fraction0["Time"].values.reshape(-1, 1), fraction0["Counts_Resistant"])
print(f"Resistant Coefficient of determination 100% (R^2): {r_sq0}")
#For b0 value, aka the intercept:
print(f"Resistant Intercept 100%: {model0.intercept_}")
#For b1 value, aka the slope
print(f"Resistant Slope 100%: {model0.coef_}")
#For Parental: 
r_sq0P = model0P.score(fraction0["Time"].values.reshape(-1, 1), fraction0["Counts_Parental"])
print(f"Parental Coefficient of determination 0% (R^2): {r_sq0P}")
#For b0 value, aka the intercept:
print(f"Parental Intercept 0%: {model0P.intercept_}")
#For b1 value, aka the slope
print(f"Parental Slope 0% : {model0P.coef_}")
# 10% Parental
r_sq10 = model10.score(fraction10["Time"].values.reshape(-1, 1), fraction10["Counts_Resistant"])
print(f"Resistant Coefficient of determination 90% (R^2): {r_sq10}")
print(f"Resistant Intercept 90%: {model10.intercept_}")
print(f"Resistant Slope 90%: {model10.coef_}")
r_sq10P = model10P.score(fraction10["Time"].values.reshape(-1, 1), fraction10["Counts_Parental"])
print(f"Parental Coefficient of determination 10% (R^2): {r_sq10P}")
print(f"Parental Intercept 10%: {model10P.intercept_}")
print(f"Parental Slope 10% : {model10P.coef_}")
# 30% Parental
r_sq30 = model30.score(fraction30["Time"].values.reshape(-1, 1), fraction30["Counts_Resistant"])
print(f"Resistant Coefficient of determination 70% (R^2): {r_sq30}")
print(f"Resistant Intercept 70%: {model30.intercept_}")
print(f"Resistant Slope 70%: {model30.coef_}")
r_sq30P = model30P.score(fraction30["Time"].values.reshape(-1, 1), fraction30["Counts_Parental"])
print(f"Parental Coefficient of determination 30% (R^2): {r_sq30P}")
print(f"Parental Intercept 30%: {model30P.intercept_}")
print(f"Parental Slope 30% : {model30P.coef_}")
# 50% Parental
r_sq50 = model50.score(fraction50["Time"].values.reshape(-1, 1), fraction50["Counts_Resistant"])
print(f"Resistant Coefficient of determination 50% (R^2): {r_sq50}")
print(f"Resistant Intercept 50%: {model50.intercept_}")
print(f"Resistant Slope 50%: {model50.coef_}")
r_sq50P = model50P.score(fraction50["Time"].values.reshape(-1, 1), fraction50["Counts_Parental"])
print(f"Parental Coefficient of determination 50% (R^2): {r_sq50P}")
print(f"Parental Intercept 50%: {model50P.intercept_}")
print(f"Parental Slope 50% : {model50P.coef_}")
# 70% Parental
r_sq70 = model70.score(fraction70["Time"].values.reshape(-1, 1), fraction70["Counts_Resistant"])
print(f"Resistant Coefficient of determination 30% (R^2): {r_sq70}")
print(f"Resistant Intercept 30%: {model70.intercept_}")
print(f"Resistant Slope 30%: {model70.coef_}")
r_sq70P = model70P.score(fraction70["Time"].values.reshape(-1, 1), fraction70["Counts_Parental"])
print(f"Parental Coefficient of determination 70% (R^2): {r_sq70P}")
print(f"Parental Intercept 70%: {model70P.intercept_}")
print(f"Parental Slope 70% : {model70P.coef_}")
# 90% Parental
r_sq90 = model90.score(fraction90["Time"].values.reshape(-1, 1), fraction90["Counts_Resistant"])
print(f"Resistant Coefficient of determination 10% (R^2): {r_sq90}")
print(f"Resistant Intercept 10%: {model90.intercept_}")
print(f"Resistant Slope 10%: {model90.coef_}")
r_sq90P = model90P.score(fraction90["Time"].values.reshape(-1, 1), fraction90["Counts_Parental"])
print(f"Parental Coefficient of determination 90% (R^2): {r_sq90P}")
print(f"Parental Intercept 90%: {model90P.intercept_}")
print(f"Parental Slope 90% : {model90P.coef_}")
# 95% Parental
r_sq95 = model95.score(fraction95["Time"].values.reshape(-1, 1), fraction95["Counts_Resistant"])
print(f"Resistant Coefficient of determination 5% (R^2): {r_sq95}")
print(f"Resistant Intercept 5%: {model95.intercept_}")
print(f"Resistant Slope 5%: {model95.coef_}")
r_sq95P = model95P.score(fraction95["Time"].values.reshape(-1, 1), fraction95["Counts_Parental"])
print(f"Parental Coefficient of determination 95% (R^2): {r_sq95P}")
print(f"Parental Intercept 95%: {model95P.intercept_}")
print(f"Parental Slope 95% : {model95P.coef_}")
# 100% Parental
r_sq1 = model1.score(fraction1["Time"].values.reshape(-1, 1), fraction1["Counts_Resistant"])
print(f"Resistant Coefficient of determination 0% (R^2): {r_sq1}")
print(f"Resistant Intercept 0%: {model1.intercept_}")
print(f"Resistant Slope 0%: {model1.coef_}")
r_sq1P = model1P.score(fraction1["Time"].values.reshape(-1, 1), fraction1["Counts_Parental"])
print(f"Parental Coefficient of determination 100% (R^2): {r_sq1P}")
print(f"Parental Intercept 100%: {model1P.intercept_}")
print(f"Parental Slope 100% : {model1P.coef_}")

#To add the linear regression equations, we need to use the numbers we have calculated above and 
#Add them to the graph. We need to also calculate the location of that we want the equation to be at 

#Function to annotate the equations within the graphs
#Within this function we are calling for a model with a number, the corresponding
#r_sq which is the R^2 value, the data which is for example Fraction0, the 
#x and y axis data which is time and counts, the label for this and the position
#In the grid, as well as the y and x position within each graph
def plot_annotate(model, r_sq, data, x, y, label, ax, ypos_factor, xpos_factor):
    # Annotate the equation
    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()
    y_pos = y_max - ypos_factor * (y_max - y_min)  # adjust position based on ypos_factor
    x_pos = x_max - xpos_factor * (x_max - x_min)
    ax.text(x_pos, y_pos, f"{label} = {model.intercept_:.2f} + {model.coef_[0]:.2f}*Time, R² = {r_sq:.2f} ", fontsize = 8)

#For fraction 0 Parental:
    #Here we use the above described function and add the appropritate values
plot_annotate(model0, r_sq0, fraction0, "Time", "Counts_Resistant", "CRES", axs[0, 0], 0.1, 0.96)
#plot_annotate(model0P, fraction0, "Time", "Counts_Parental", "CPAR", axs[0, 0], 0.2)
#For Fraction 10% Parental
plot_annotate(model10, r_sq10, fraction10, "Time", "Counts_Resistant", "CRES", axs[1, 0], 0.1, 0.96)
plot_annotate(model10P, r_sq10P, fraction10, "Time", "Counts_Parental", "CPAR", axs[1, 0], 0.2, 0.96)
#For Fraction 30% Parental
plot_annotate(model30, r_sq30, fraction30, "Time", "Counts_Resistant", "CRES", axs[2, 0], 0.8, 0.59)
plot_annotate(model30P, r_sq30P, fraction30, "Time", "Counts_Parental", "CPAR", axs[2, 0], 0.9, 0.59)
#For Fraction 50% Parental
plot_annotate(model50, r_sq50, fraction50, "Time", "Counts_Resistant", "CRES", axs[3, 0], 0.8, 0.59)
plot_annotate(model50P, r_sq50P, fraction50, "Time", "Counts_Parental", "CPAR", axs[3, 0], 0.9, 0.59)
#For Fraction 70% Parental
plot_annotate(model70, r_sq70, fraction70, "Time", "Counts_Resistant", "CRES", axs[0, 1], 0.8, 0.59)
plot_annotate(model70P, r_sq70P, fraction70, "Time", "Counts_Parental", "CPAR", axs[0, 1], 0.9, 0.59)
#For Fraction 90% Parental
plot_annotate(model90, r_sq90, fraction90, "Time", "Counts_Resistant", "CRES", axs[1, 1], 0.8, 0.59)
plot_annotate(model90P, r_sq90P, fraction90, "Time", "Counts_Parental", "CPAR", axs[1, 1], 0.9, 0.59)
#For Fraction 95% Parental
plot_annotate(model95, r_sq95, fraction95, "Time", "Counts_Resistant", "CRES", axs[2, 1], 0.8, 0.59)
plot_annotate(model95P, r_sq95P, fraction95, "Time", "Counts_Parental", "CPAR", axs[2, 1], 0.9, 0.59)
#For Fraction 100% Parental
#plot_annotate(model1, r_sq1, fraction1, "Time", "Counts_Resistant", "CRES", axs[3, 1], 0.8, 0.59)
plot_annotate(model1P, r_sq1P, fraction1, "Time", "Counts_Parental", "CPAR", axs[3, 1], 0.9, 0.59)

plt.subplots_adjust(hspace=0.3)
plt.suptitle("Counts at different parental fractions")
plt.show()

#Now i want to display the slope values in a separate table which is plotted
#instead of just outputted
# We need to create a list of lists by using the following:
# Here we need to give a label using"" and the actual value which is calculated earlier
data = [["Parental 0%", model0P.coef_],
        ["Resistant 100%", model0.coef_],
        ["Parental 10%", model10P.coef_],
        ["Resistant 90%", model10.coef_],
        ["Parental 30%", model30P.coef_],
        ["Resistant 70%", model30.coef_], 
        ["Parental 50%", model50P.coef_],
        ["Resistant 50%", model50.coef_],
        ["Parental 70%", model70P.coef_],
        ["Resistant 30%", model70.coef_],
        ["Parental 90%", model90P.coef_],
        ["Resistant 10%", model90.coef_],
        ["Parental 95%", model95P.coef_],
        ["Resistant 5%", model95.coef_], 
        ["Parental 100%", model1P.coef_],
        ["Resistant 0%", model1.coef_]]
       


# Create a figure and axes
fig, ax = plt.subplots()

# Hide axes
ax.axis('off')

# Create a table and add it to the axes
table = plt.table(cellText=data, 
                  colLabels=["Fraction at start", "Slope"], 
                  cellLoc = 'center', 
                  loc='center')

# Scale the table
table.scale(1, 1.5)

plt.show()
#Also want to make one for the R^2 of each linear regression line
data = [["Parental 0%", r_sq0P],
        ["Resistant 100%", r_sq0],
        ["Parental 10%", r_sq10P],
        ["Resistant 90%", r_sq10],
        ["Parental 30%", r_sq30P],
        ["Resistant 70%", r_sq30], 
        ["Parental 50%", r_sq50P],
        ["Resistant 50%", r_sq50],
        ["Parental 70%", r_sq70P],
        ["Resistant 30%", r_sq70],
        ["Parental 90%", r_sq90P],
        ["Resistant 10%", r_sq90],
        ["Parental 95%", r_sq95P],
        ["Resistant 5%", r_sq95], 
        ["Parental 100%", r_sq1P],
        ["Resistant 0%", r_sq1]]
       


# Create a figure and axes
fig, ax = plt.subplots()

# Hide axes
ax.axis('off')

# Create a table and add it to the axes
table = plt.table(cellText=data, 
                  colLabels=["fraction at start", "R² Value"], 
                  cellLoc = 'center', 
                  loc='center')

# Scale the table
table.scale(1, 1.5)