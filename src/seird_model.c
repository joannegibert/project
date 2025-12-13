#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*****************************************************************
Workflow:
1. Read multiple sets of model parameters (beta, gamma, mu, N, sigma)
   from an input CSV file, where each row represents a different scenario.
2. For each scenario:
   - Initialize the population compartments (S, E, I, R, D).
   - Simulate the epidemic over a number of time steps using the Euler method.
   - Record the values of S, E, I, R, D at each step.
   - Save the results to a separate output CSV file.
3. Repeat for all scenarios.

Outcome:
Each scenario produces a CSV file showing how the epidemic evolves over time,
allowing for analysis and comparison in Matlab.
******************************************************************/

//   -------------- MODEL PARAMETERS --------------
// This structure groups all parameters needed by the SEIRD model.
// Instead of passing many variables to functions, we pass a single
// struct (p) containing every parameter.

struct Params {
    double beta;   // transmission rate [day⁻¹]
    double sigma;  // incubation rate [day⁻¹]
    double gamma;  // recovery rate [day⁻¹]
    double mu;     // mortality rate [day⁻¹]
    double N;      // total population
    double dt;     // time step [days]
    int steps;     // total simulation steps
};

//   -------------- SEIRD DIFFERENTIAL EQUATIONS + EULER UPDATE --------------
// Function to perform one Euler time-step: updates S, E, I, R, D from time t to t + dt
// Inputs:
//   - S, E, I, R, D : pointers to the current state variables of the SEIRD model
//     (Susceptible, Exposed, Infected, Recovered, Dead). Their values at time t are read and updated in place.
//   - p : Params structure containing all model parameters
//     (beta, sigma, gamma, mu, dt, N)
//
// Output:
//   - Updates *S, *E, *I, *R, *D using one Euler integration step.
//     No return value: the function directly modifies the variables via pointers.
//  --------------------------------------------------------------------------

void seird_step(double *S, double *E, double *I, double *R, double *D, struct Params p) {
    double dS = -p.beta * (*S) * (*I) / p.N;                  // dS/dt = −β SI/N
    double dE =  p.beta * (*S) * (*I) / p.N - p.sigma * (*E); // dE/dt = β SI/N −σE
    double dI =  p.sigma * (*E) - (p.gamma + p.mu) * (*I);    // dI/dt = σE −(γ+μ)I
    double dR =  p.gamma * (*I);                              // dR/dt = γI
    double dD =  p.mu * (*I);                                 // dD/dt = μI

    // Euler update over time: X(t + dt) = X(t) + (dX/dt) * dt
    // here dX variables represent time derivatives dX/dt evaluated at the current state
    *S += dS * p.dt;
    *E += dE * p.dt;
    *I += dI * p.dt;
    *R += dR * p.dt;
    *D += dD * p.dt;
}

//   -------------- RUN SEIRD MODEL + CREATE CSV --------------
// Function to run the SEIRD simulation using the given parameters and store the results to a CSV file

// Inputs:
//   - p: Params structure containing all model parameters (N, beta, gamma,
//     sigma, mu, dt, steps)
//   - filename: name of the CSV file to create and write results into
//
// Output:
//   - creates a CSV file containing the time evolution of the SEIRD model
//     (Columns: time, S, E, I, R, D)
//   - prints a confirmation message when the file is successfully written and saved
// -------------------------------------------------------------

void run_seird_model(struct Params p, const char *filename) {

     // Open the output CSV file
    FILE *f = fopen(filename, "w"); // w = write mode
    if (!f) {
        printf("Error: could not open %s\n", filename);
        exit(1);
    }
    
    // ----- INITIAL CONDITIONS -----
    double I0 = 1; // Start with 1 infected individual
    double S = p.N - I0;
    double E = 0;
    double I = I0;
    double R = 0;
    double D = 0;

    // Write CSV header
    fprintf(f, "time,S,E,I,R,D\n");

    // ----- TIME LOOP -----
    // For each time step, save the current values then update them
    for (int t = 0; t < p.steps; t++) {

        double time = t * p.dt; // Convert step number into time in days
        
        // Write current state to CSV file
        fprintf(f, "%.2f,%.4f,%.4f,%.4f,%.4f,%.4f\n", time, S, E, I, R, D);
        
        // Update S,E,I,R,D  
        // Pass the addresses of S, E, I, R, D so seird_step can modify their values directly
        seird_step(&S, &E, &I, &R, &D, p);
    }

    fclose(f);
    printf("File %s created.\n", filename);
}

//   -------------- MAIN FUNCTION --------------
// Reads parameters from "parameters.csv", runs one simulation
// for each row of the file, and generates one output CSV per scenario
int main() {

    // Open the input CSV file containing the parameters
    FILE *csv = fopen("../data/parameters.csv", "r"); // r = read mode
    if (!csv) {
        printf("Error: could not open parameters.csv\n");
        return 1;
    }

    struct Params p;
    p.dt = 0.1;         // time step (days)
    p.steps = 20000;    // number of steps (2000 days)

    char line[256];     // Store the content of the CSV file
    int line_count = 0; // keeps track of which line (e.g. which scenario) we are currently reading from the CSV file

    // Read the file line by line
    while (fgets(line, sizeof(line), csv)) {
        // Skip the header (line 0)
        if (line_count == 0) { 
            line_count++; 
            continue; 
        }

        // Extract beta, gamma, mu, N, sigma from the line (as double)
        sscanf(line, "%lf,%lf,%lf,%lf,%lf",
                   &p.beta, &p.gamma, &p.mu, &p.N, &p.sigma)

            // Create a filename like "scenario_1.csv", "scenario_2.csv", ...
            char filename[64];
            sprintf(filename, "../results/scenario_%d.csv", line_count);

            // Run SEIRD model with parameters from this row
            run_seird_model(p, filename);

        line_count++;
    }

    fclose(csv);
    printf("All simulations completed (%d scenarios)\n", line_count - 1); // line_count - header 
    return 0;
}
