import csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

avg_wavelength = []
cross_section = []

with open('jpl_cm2.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        
        waveln_limits = row["lambda"].split("–")
        avg_wavelength.append( .5*(float(waveln_limits[0]) + float(waveln_limits[1])))
        cross_section.append(row["k294"])

avg_wavelength = np.array(avg_wavelength).reshape((-1, 1))
cross_section = np.array(cross_section)
xs_model = LinearRegression().fit(avg_wavelength, cross_section)


# Create quantum_yield dictionary
wavelength_qy = {}
with open('quantum_yield.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        wavelength_qy[float(row["wavelength"])] = float(row["k298"])

# calculate the constant for each zenith angle
angles = [0, 10 ,20,30,40,50,60,70,78,86]
angles_ks = {}
for angle in angles:
    angles_ks[angle] = 0

with open('actinic flux.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        # calculate average wavelength

        wavelengths = row["Wavelength"].split("-")
        avg_wavelength = .5* (float(wavelengths[0]) + float(wavelengths[1]))
    
        xsec = xs_model.intercept_ + xs_model.coef_*avg_wavelength
        
        qy = 1
        if avg_wavelength > 398 and avg_wavelength <= 422:
            qy = wavelength_qy[int(avg_wavelength) + 1]
            
        elif avg_wavelength > 422:
            qy = 0
        
        for angle in angles:
            #actinic flux already multiplied by delta * lambda
            # Multiply by 60 to convert to min-1
            k_tmp = (10**int(row["Exponent"])) * float(row[str(angle)]) *(10**-20) * xsec[0] *  qy * 90

            angles_ks[angle] += k_tmp
    
    # Calculate polynomial fit
    ks = []
    for angle in angles:
        print(angles_ks[angle])
        ks.append(angles_ks[angle])

    # Plot rate coefficients vs. cosine of zenith angle
    cos_angles = [np.cos(np.radians(a)) for a in angles]
    plt.figure()
    plt.plot(cos_angles, ks, "or")
    plt.xlabel("cosine of zenith angle")
    plt.ylabel("kNO2 (min-1)")
    plt.title("kNO2 (min -1) vs cosine of zenith angle")
    plt.show()
    plt.figure()

    coefficients = np.polyfit(np.array(cos_angles), np.array(ks), 4)
    # print(coefficients)
    p = np.poly1d(coefficients)


    predicted = []
    t = []
    actual = []
    with open('experimental_data.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        
        for row in reader:

            
            my_angle_cos = np.cos(np.radians(float(row["angle"])))
            expected_k = p(my_angle_cos)
            
            predicted.append(expected_k)
            t.append(int(row["time"]))
            if float(row["kNO2"]) != 0:
                actual.append(float(row["kNO2"]))
            else:
                actual.append(np.nan)

    # Create the plot
    plt.plot(t, actual, "or", label = "Experimental")
    plt.plot(t, predicted, label = "Predicted")
    plt.xlabel("Time (PST)")
    plt.ylabel("kNO2 (min-1)")
    ax = plt.gca()
    ax.legend(['Experimental', 'Predicted'])
    plt.title("Experimental vs. Predicted values for NO2 photolysis in El Monte, CA")
    plt.show()




