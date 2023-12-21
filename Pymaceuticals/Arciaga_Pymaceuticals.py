# Arciaga - Pymaceuticals Inc. 
---

### Analysis 
### 1. Summary Statistics: 
- The summary statistics table provides key insights into the distribution of tumor volumes for each drug regimen. Capomulin and Ramicane show the lowest mean and median tumor volumes, which suggests the potential effectiveness in reducing tumor sizes compared to other drug regimens in this dataset. 
- This table also provides the standard error of the mean (SEM), which gives an indication of result reliability. Capomulin and Ramicane have the lowest SEM values, suggesting more precise estimates. 
### 2. Tumor Volume Distribution: 
- The boxplot illustrated the distribution of final tumor volumes across the selected drug regimens: Capomulin, Ramicane, Infubinol, and Ceftamin. Capomulin and Ramicane show lower median final tumor volumes and smaller interquartile ranges, indicating a more consistent and potentially effective treatment compared to Infubinol and Ceftamin.  
- Outlier analysis identifies potential outliers in the Infubinol regimen, suggesting some variability in treatment response within this group.  
### 3. Correlation and Regression Analysis
- The scatter plot and linear regression analysis for mouse weight vs. average tumor volume under the Capomulin regimen displays a positive correlation. As mouse weight increases, the average tumor volume tends to increase as well. The regression line further supports this trend, indicating a potential positive relationship between weight and treatment effectiveness.  
- The correlation coefficient of approximately 0.84 strengthens the evidence of a positive correlation. This implies there is a direct relationship between mouse weight and the observed tumor volume in response to the Capomulin regimen. 

# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single DataFrame
mouse_data_complete = pd.merge(study_results, mouse_metadata, how="left", on=["Mouse ID", "Mouse ID"])

# Display the data table for preview
mouse_data_complete.head()

# Checking the number of mice.
mice_count = mouse_data_complete["Mouse ID"].nunique()
mice_count

# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mice = mouse_data_complete[mouse_data_complete.duplicated(subset=["Mouse ID", "Timepoint"])]["Mouse ID"].unique()
duplicate_mice

# Optional: Get all the data for the duplicate mouse ID. 
mouse_data_complete.loc[mouse_data_complete["Mouse ID"]=="g989"]

# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_mouse_data_complete = mouse_data_complete.loc[mouse_data_complete["Mouse ID"]!="g989"]
clean_mouse_data_complete

# Checking the number of mice in the clean DataFrame.
clean_mouse_data_complete["Mouse ID"].nunique()

# Arciaga - Summary Statistics

# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
means = clean_mouse_data_complete.groupby("Drug Regimen")["Tumor Volume (mm3)"].mean()
medians = clean_mouse_data_complete.groupby("Drug Regimen")["Tumor Volume (mm3)"].median()
variances = clean_mouse_data_complete.groupby("Drug Regimen")["Tumor Volume (mm3)"].var()
stds = clean_mouse_data_complete.groupby("Drug Regimen")["Tumor Volume (mm3)"].std()
sems = clean_mouse_data_complete.groupby("Drug Regimen")["Tumor Volume (mm3)"].sem()

summary_statistics_table = pd.DataFrame({
    "Mean Tumor Volume":means,
    "Median Tumor Volume":medians,
    "Tumor Volume Variance":variances,
    "Tumor Volume Std. Dev.":stds,
    "Tumor Volume Std. Err.":sems
})
summary_statistics_table

# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
summary_statistics_adv = clean_mouse_data_complete.groupby("Drug Regimen").agg({"Tumor Volume (mm3)":["mean","median","var","std","sem"]})
summary_statistics_adv

# Arciaga - Bar and Pie Charts

# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
drug_regimen_count = clean_mouse_data_complete["Drug Regimen"].value_counts()

drug_regimen_count.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("# of Observed Mouse Timepoints")
plt.show()

# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
drug_regimen_count = clean_mouse_data_complete["Drug Regimen"].value_counts()

plt.bar(drug_regimen_count.index.values, drug_regimen_count.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("# of Observed Mouse Timepoints")
plt.show()

# Generate a pie plot showing the distribution of female versus male mice using Pandas
sex_count = clean_mouse_data_complete["Sex"].value_counts()

sex_count.plot(kind="pie",autopct='%1.1f%%')
plt.ylabel("Sex")
plt.show()

# Generate a pie plot showing the distribution of female versus male mice using pyplot
sex_count = clean_mouse_data_complete["Sex"].value_counts()

plt.pie(sex_count, labels=sex_count.index, autopct="%1.1f%%",)
plt.title("Distribution of Female vs. Male Mice")
plt.ylabel("Sex")
plt.show()

# Arciaga - Quartiles, Outliers and Boxplots

# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
max_timepoints = clean_mouse_data_complete.groupby("Mouse ID")["Timepoint"].max()
final_tumor_volume = max_timepoints.reset_index()

#Select data for the four specified regims
selected_regimens = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
final_tumor_volume_grouped = final_tumor_volume.merge(clean_mouse_data_complete, on = ["Mouse ID", "Timepoint"], how = "left")
final_tumor_volume_grouped

# Put treatments into a list for for loop (and later for plot labels)
treatments = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_data = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for treatment in treatments:
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    tumor_vol = final_tumor_volume_grouped.loc[final_tumor_volume_grouped["Drug Regimen"] == treatment, "Tumor Volume (mm3)"]
    
    # add subset 
    tumor_vol_data.append(tumor_vol)
    
    # Calculate the IQR for each treatment regimen
    quartiles = tumor_vol.quantile([0.25, 0.5, 0.75])
    lower_q = quartiles[0.25]
    upper_q = quartiles[0.75]
    iqr = upper_q - lower_q

    # Determine upper and lower bounds to identify potential outliers
    lower_bound = lower_q - 1.5 * iqr
    upper_bound = upper_q + 1.5 * iqr
    
    # Determine outliers using upper and lower bounds
    potential_outliers = tumor_vol[(tumor_vol < lower_bound) | (tumor_vol > upper_bound)]
    
    # Print results
    print(f"{treatment}'s potential outliers: {potential_outliers}")
    
# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
plt.figure(figsize=(5, 5))
plt.boxplot(tumor_vol_data, labels=treatments, sym='ro')
plt.title("Tumor Volume Across Selected Regimens with Potential Outliers")
plt.xlabel("Drug Regimen")
plt.ylabel("Final Tumor Volume (mm3)")
plt.show()

# Arciaga - Line and Scatter Plots

# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
# Choose the first mouse in the dataset that was treated with Capomulin
selected_mouse_id = clean_mouse_data_complete.loc[clean_mouse_data_complete["Drug Regimen"] == "Capomulin", "Mouse ID"].iloc[0]

# Filter the data for the selected mouse and Capomulin regimen
selected_mouse_data = clean_mouse_data_complete.loc[(clean_mouse_data_complete["Mouse ID"] == selected_mouse_id) & (clean_mouse_data_complete["Drug Regimen"] == "Capomulin")]

plt.figure(figsize=(5, 5))
plt.plot(selected_mouse_data["Timepoint"], selected_mouse_data["Tumor Volume (mm3)"], linestyle='-',)
plt.title(f"Capomulin treatment of mouse {selected_mouse_id}")
plt.xlabel("Timepoint")
plt.ylabel("Tumor Volume (mm3)")
plt.show()

# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
# Calculate the average tumor volume for each mouse in the Capomulin regimen
avg_tumor_volume_capomulin = clean_mouse_data_complete[clean_mouse_data_complete["Drug Regimen"] == "Capomulin"].groupby("Mouse ID")["Tumor Volume (mm3)"].mean()

# Merge the average tumor volume data with the mouse metadata to get the weight information
capomulin_data = pd.merge(avg_tumor_volume_capomulin, mouse_metadata, on="Mouse ID")
print(capomulin_data)

# Create a scatter plot
plt.figure(figsize=(5, 5))
plt.scatter(capomulin_data["Weight (g)"], capomulin_data["Tumor Volume (mm3)"], marker="o", alpha=0.75)
plt.title("Mouse Weight vs. Average Tumor Volume for Capomulin Regimen")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.show()

# Arciaga - Correlation and Regression

# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen
correlation_coefficient = st.pearsonr(capomulin_data["Weight (g)"], capomulin_data["Tumor Volume (mm3)"])[0]
print(f"The correlation coefficient between mouse weight and average tumor volume is {correlation_coefficient:.2f}")

# Calculate the linear regression model for mouse weight and average tumor volume for Capomulin regimen
(slope, intercept, rvalue, pvalue, stderr) = st.linregress(capomulin_data["Weight (g)"], capomulin_data["Tumor Volume (mm3)"])

# Create the regression equation
regression_equation = f"y = {slope:.2f}x + {intercept:.2f}"

# Create the scatter plot
plt.figure(figsize=(5, 5))
plt.scatter(capomulin_data["Weight (g)"], capomulin_data["Tumor Volume (mm3)"], marker="o", alpha=0.75, label="Data Points")

# Plot the regression line
regress_values = slope * capomulin_data["Weight (g)"] + intercept
plt.plot(capomulin_data["Weight (g)"], regress_values, "r-", label=regression_equation)

# Add labels and a legend
plt.title("Mouse Weight vs. Avg Tumor Vol. for Capomulin Regimen")
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.legend()

# Show the plot
plt.show()
