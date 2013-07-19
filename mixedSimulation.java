package antibioticResistance;

import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: ellscampbell
 * Date: 7/16/13
 * Time: 4:36 PM
 * To change this template use File | Settings | File Templates.
 */
public class mixedSimulation {

    // variable reference and storage
    private static final int micA = 0;
    private static final int micB = 1;

    private static final int uninfected = 0;

    private static final int susceptible_drugA = 1;
    private static final int susceptible_drugB = 2;

    private static final int resistantA_drugA = 3;
    private static final int resistantA_drugB = 4;

    private static final int resistantB_drugA = 5;
    private static final int resistantB_drugB = 6;

    private static final int multiResistant_drugA = 7;
    private static final int multiResistant_drugB = 8;

    double[] pharmacoGrowth = new double[9];
    double[] patients = new double[9];
    private int deltaTime = 1;
    private double PsiMax, PsiMin, mu, immuneClearance, rateOfInfection, rateOfSuperInfection, cost;
    private int dosageA, dosageB;
    private int period, burn, timeMax;

    ArrayList<double[]> frequencies = new ArrayList<double[]>();

    public static void main(String[] args) {
        mixedSimulation mix = new mixedSimulation();
        mix.run();
    }

    public void run() {
        initGrowth();
        for (int timeStep = 1; timeStep < timeMax; timeStep++) {

            if (timeStep%period == 0) {
                switchPopulations();
            }

            fourthOrderApproximation();
            frequencies.add(patients);

        }

        findMeans();
    }


    private void switchPopulations() {

        // switch S



        double[] holderS = new double[2];
        holderS[0] = patients[susceptible_drugA];
        holderS[1] = patients[susceptible_drugB];

        patients[susceptible_drugA] = 0;
        patients[susceptible_drugB] = 0;

        patients[susceptible_drugA] = holderS[1];
        patients[susceptible_drugB] = holderS[0];



        // switch R_drugA
        double[] holderRA = new double[2];
        holderRA[0] = patients[resistantA_drugA];
        holderRA[1] = patients[resistantA_drugB];

        patients[resistantA_drugA] = 0;
        patients[resistantA_drugB] = 0;

        patients[resistantA_drugA] = holderRA[1];
        patients[resistantA_drugB] = holderRA[0];

        // switch R_drugB
        double[] holderRB = new double[2];
        holderRB[0] = patients[resistantB_drugA];
        holderRB[1] = patients[resistantB_drugB];

        patients[resistantB_drugA] = 0;
        patients[resistantB_drugB] = 0;

        patients[resistantB_drugA] = holderRB[1];
        patients[resistantB_drugB] = holderRB[0];

        // switch R_multiDrug
        double[] holderRM = new double[2];
        holderRM[0] = patients[multiResistant_drugA];
        holderRM[1] = patients[multiResistant_drugB];

        patients[multiResistant_drugA] = 0;
        patients[multiResistant_drugB] = 0;

        patients[multiResistant_drugA] = holderRM[1];
        patients[multiResistant_drugB] = holderRM[0];

        initGrowth();

    }

    public void initConstants() {
        // bergstrom constants
        cost = SimulationSettings.getInstance().getCost();
        immuneClearance = SimulationSettings.getInstance().getImmuneClearance();
        rateOfInfection = SimulationSettings.getInstance().getRateOfInfection();
        rateOfSuperInfection = SimulationSettings.getInstance().getRateOfSuperInfection();
        period = SimulationSettings.getInstance().getPeriod();

        // pharmacodynamic constants
        PsiMax = SimulationSettings.getInstance().getMaximumBacterialGrowthRate();
        PsiMin = SimulationSettings.getInstance().getMinimumBacterialGrowthRate_WithMaxDrug();

        // our constants
        mu = SimulationSettings.getInstance().getMutationRate();

        burn = period * 10000;
        timeMax = period * 100000;

    }

    private void fourthOrderApproximation() {
        ArrayList<double[]> fourDeltas = new ArrayList<double[]>() ;

        for (int i = 0; i < 3; i++) {
            double[] deltaPatients = new double[9];
            deltaPatients[uninfected] = deltaUninfected();
            deltaPatients[susceptible_drugA] = deltaSusceptible_drugA();
            deltaPatients[susceptible_drugB] = deltaSusceptible_drugB();

            deltaPatients[resistantA_drugA] = deltaResistantA_drugA();
            deltaPatients[resistantA_drugB] = deltaResistantA_drugB();

            deltaPatients[resistantB_drugA] = deltaResistantB_drugA();
            deltaPatients[resistantB_drugB] = deltaResistantB_drugB();

            deltaPatients[multiResistant_drugA] = deltaMultiResistant_drugA();
            deltaPatients[multiResistant_drugB] = deltaMultiResistant_drugB();

            patients[uninfected]     += 0.5 * deltaTime * deltaPatients[uninfected];

            patients[susceptible_drugA]    += 0.5 * deltaTime * deltaPatients[susceptible_drugA];
            patients[susceptible_drugB]    += 0.5 * deltaTime * deltaPatients[susceptible_drugB];

            patients[resistantA_drugA]     += 0.5 * deltaTime * deltaPatients[resistantA_drugA];
            patients[resistantA_drugB]     += 0.5 * deltaTime * deltaPatients[resistantA_drugB];

            patients[resistantB_drugA]     += 0.5 * deltaTime * deltaPatients[resistantB_drugA];
            patients[resistantB_drugB]     += 0.5 * deltaTime * deltaPatients[resistantB_drugB];

            patients[multiResistant_drugA] += 0.5 * deltaTime * deltaPatients[multiResistant_drugA];
            patients[multiResistant_drugB] += 0.5 * deltaTime * deltaPatients[multiResistant_drugB];


            fourDeltas.add(deltaPatients);
        }

        double[] deltaPatients = new double[9];

        deltaPatients[uninfected] = deltaTime * deltaUninfected();
        deltaPatients[susceptible_drugA] = deltaTime * deltaSusceptible_drugA();
        deltaPatients[susceptible_drugB] = deltaTime * deltaSusceptible_drugB();

        deltaPatients[resistantA_drugA] = deltaTime * deltaResistantA_drugA();
        deltaPatients[resistantA_drugB] = deltaTime * deltaResistantA_drugB();

        deltaPatients[resistantB_drugA] = deltaTime * deltaResistantB_drugA();
        deltaPatients[resistantB_drugB] = deltaTime * deltaResistantB_drugB();

        deltaPatients[multiResistant_drugA] = deltaTime * deltaMultiResistant_drugA();
        deltaPatients[multiResistant_drugB] = deltaTime * deltaMultiResistant_drugB();

        fourDeltas.add(deltaPatients);



        patients[uninfected]     += (fourDeltas.get(0)[uninfected] + (2.0 * fourDeltas.get(1)[uninfected]) + (2.0 * fourDeltas.get(2)[uninfected]) + fourDeltas.get(3)[uninfected]) * ((1.0/6.0) * deltaTime);

        patients[susceptible_drugA]    += (fourDeltas.get(0)[susceptible_drugA] + (2.0 * fourDeltas.get(1)[susceptible_drugA]) + (2.0 * fourDeltas.get(2)[susceptible_drugA]) + fourDeltas.get(3)[susceptible_drugA]) * ((1.0/6.0) * deltaTime);
        patients[susceptible_drugB]    += (fourDeltas.get(0)[susceptible_drugB] + (2.0 * fourDeltas.get(1)[susceptible_drugB]) + (2.0 * fourDeltas.get(2)[susceptible_drugB]) + fourDeltas.get(3)[susceptible_drugB]) * ((1.0/6.0) * deltaTime);

        patients[resistantA_drugA]     += (fourDeltas.get(0)[resistantA_drugA] + (2.0 * fourDeltas.get(1)[resistantA_drugA]) + (2.0 * fourDeltas.get(2)[resistantA_drugA]) + fourDeltas.get(3)[resistantA_drugA]) * ((1.0/6.0) * deltaTime);
        patients[resistantA_drugB]     += (fourDeltas.get(0)[resistantA_drugB] + (2.0 * fourDeltas.get(1)[resistantA_drugB]) + (2.0 * fourDeltas.get(2)[resistantA_drugB]) + fourDeltas.get(3)[resistantA_drugB]) * ((1.0/6.0) * deltaTime);

        patients[resistantB_drugA]     += (fourDeltas.get(0)[resistantB_drugA] + (2.0 * fourDeltas.get(1)[resistantB_drugA]) + (2.0 * fourDeltas.get(2)[resistantB_drugA]) + fourDeltas.get(3)[resistantB_drugA]) * ((1.0/6.0) * deltaTime);
        patients[resistantB_drugB]     += (fourDeltas.get(0)[resistantB_drugB] + (2.0 * fourDeltas.get(1)[resistantB_drugB]) + (2.0 * fourDeltas.get(2)[resistantB_drugB]) + fourDeltas.get(3)[resistantB_drugB]) * ((1.0/6.0) * deltaTime);

        patients[multiResistant_drugA] += (fourDeltas.get(0)[multiResistant_drugA] + (2.0 * fourDeltas.get(1)[multiResistant_drugA]) + (2.0 * fourDeltas.get(2)[multiResistant_drugA]) + fourDeltas.get(3)[multiResistant_drugA]) * ((1.0/6.0) * deltaTime);
        patients[multiResistant_drugB] += (fourDeltas.get(0)[multiResistant_drugB] + (2.0 * fourDeltas.get(1)[multiResistant_drugB]) + (2.0 * fourDeltas.get(2)[multiResistant_drugB]) + fourDeltas.get(3)[multiResistant_drugB]) * ((1.0/6.0) * deltaTime);



    }

    public void initFrequencies() {
        this.patients[uninfected] = SimulationSettings.getInstance().getHealthyPatients();

        this.patients[susceptible_drugA] = SimulationSettings.getInstance().getSusceptible_InfectedPatients()/ 2.0 ;
        this.patients[susceptible_drugB] = SimulationSettings.getInstance().getSusceptible_InfectedPatients()/ 2.0 ;

        this.patients[resistantA_drugA] = SimulationSettings.getInstance().getResistantA_InfectedPatients()/ 2.0 ;
        this.patients[resistantA_drugB] = SimulationSettings.getInstance().getResistantA_InfectedPatients()/ 2.0 ;

        this.patients[resistantB_drugA] = SimulationSettings.getInstance().getResistantB_InfectedPatients()/ 2.0 ;
        this.patients[resistantB_drugB] = SimulationSettings.getInstance().getResistantB_InfectedPatients()/ 2.0 ;

        this.patients[multiResistant_drugA] = SimulationSettings.getInstance().getMultiResistant_InfectedPatients()/ 2.0;
        this.patients[multiResistant_drugB] = SimulationSettings.getInstance().getMultiResistant_InfectedPatients()/ 2.0;
    }


    public void initGrowth() {
        this.pharmacoGrowth[uninfected]           = 0;

        this.pharmacoGrowth[susceptible_drugA]    = calculateNetGrowth(susceptible_drugA);
        this.pharmacoGrowth[susceptible_drugB]    = calculateNetGrowth(susceptible_drugB);

        this.pharmacoGrowth[resistantA_drugA]     = calculateNetGrowth(resistantA_drugA);
        this.pharmacoGrowth[resistantA_drugB]     = calculateNetGrowth(resistantA_drugB);

        this.pharmacoGrowth[resistantB_drugA]     = calculateNetGrowth(resistantB_drugA);
        this.pharmacoGrowth[resistantB_drugB]     = calculateNetGrowth(resistantB_drugB);

        this.pharmacoGrowth[multiResistant_drugA] = calculateNetGrowth(multiResistant_drugA);
        this.pharmacoGrowth[multiResistant_drugB] = calculateNetGrowth(multiResistant_drugB);

    }


    private double calculateNetGrowth(int strain) {
        double[] MICs = selectMICs(strain);
        double micA = MICs[this.micA];
        double micB = MICs[this.micB];
        double netGrowth = 0;

        double effectA = ((PsiMax - PsiMin) * ((dosageA/ micA ))) / (((dosageA/ micA )) - (PsiMin / PsiMax));
        double effectB = ((PsiMax - PsiMin) * ((dosageB/ micB ))) / (((dosageB/ micB )) - (PsiMin / PsiMax));

        netGrowth = PsiMax - effectA - effectB;

        return netGrowth;
    }

    private double[] selectMICs(int strain) {
        double[] MICs = new double[2];
        double micA = 0;
        double micB = 0;

        if (strain == susceptible_drugA) {

            micA = SimulationSettings.getInstance().getMicSusceptible_drugA();
            micB = SimulationSettings.getInstance().getMicSusceptible_drugB();
            SimulationSettings.getInstance().setSingleDosageA(240);
            SimulationSettings.getInstance().setSingleDosageB(0);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();

        }

        if (strain == susceptible_drugB) {
            micA = SimulationSettings.getInstance().getMicSusceptible_drugA();
            micB = SimulationSettings.getInstance().getMicSusceptible_drugB();
            SimulationSettings.getInstance().setSingleDosageA(0);
            SimulationSettings.getInstance().setSingleDosageB(240);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        if (strain == resistantA_drugA) {
            micA = SimulationSettings.getInstance().getMicResistantA_drugA();
            micB = SimulationSettings.getInstance().getMicResistantA_drugB();
            SimulationSettings.getInstance().setSingleDosageA(240);
            SimulationSettings.getInstance().setSingleDosageB(0);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        if (strain == resistantA_drugB) {
            micA = SimulationSettings.getInstance().getMicResistantA_drugA();
            micB = SimulationSettings.getInstance().getMicResistantA_drugB();
            SimulationSettings.getInstance().setSingleDosageA(0);
            SimulationSettings.getInstance().setSingleDosageB(240);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        if (strain == resistantB_drugA) {
            micA = SimulationSettings.getInstance().getMicResistantB_drugA();
            micB = SimulationSettings.getInstance().getMicResistantB_drugB();
            SimulationSettings.getInstance().setSingleDosageA(240);
            SimulationSettings.getInstance().setSingleDosageB(0);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        if (strain == resistantB_drugB) {
            micA = SimulationSettings.getInstance().getMicResistantB_drugA();
            micB = SimulationSettings.getInstance().getMicResistantB_drugB();
            SimulationSettings.getInstance().setSingleDosageA(0);
            SimulationSettings.getInstance().setSingleDosageB(240);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        if (strain == multiResistant_drugA) {
            micA = SimulationSettings.getInstance().getMicMultiResistant_drugA();
            micB = SimulationSettings.getInstance().getMicMultiResistant_drugB();
            SimulationSettings.getInstance().setSingleDosageA(240);
            SimulationSettings.getInstance().setSingleDosageB(0);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        if (strain == multiResistant_drugB) {
            micA = SimulationSettings.getInstance().getMicMultiResistant_drugA();
            micB = SimulationSettings.getInstance().getMicMultiResistant_drugB();
            SimulationSettings.getInstance().setSingleDosageA(0);
            SimulationSettings.getInstance().setSingleDosageB(240);
            dosageA = SimulationSettings.getInstance().getSingleDosageA();
            dosageB = SimulationSettings.getInstance().getSingleDosageB();
        }

        MICs[this.micA] = micA;
        MICs[this.micB] = micB;
        return MICs;
    }


    private double deltaSusceptible_drugA() {
        double deltaS = 0;
        double mutation, clearance, infection, superInfection;


        mutation =  (mu * patients[resistantA_drugA]) + (mu * patients[resistantB_drugA]) + (mu * patients[multiResistant_drugA]) - (3 * mu * patients[susceptible_drugA]);
        clearance =  (pharmacoGrowth[susceptible_drugA] - immuneClearance) * patients[susceptible_drugA];
        infection =  ((patients[susceptible_drugA] + patients[susceptible_drugB])/2.0) * patients[uninfected] * rateOfInfection;
        superInfection = ((patients[susceptible_drugA] + patients[susceptible_drugB]) * rateOfSuperInfection * rateOfInfection * cost) * (patients[resistantA_drugA] + patients[resistantB_drugA] + patients[multiResistant_drugA]);

        deltaS = mutation + clearance + infection + superInfection ;

        return deltaS;
    }

    private double deltaSusceptible_drugB() {
        double deltaS = 0;
        double mutation, clearance, infection, superInfection;


        mutation =  (mu * patients[resistantA_drugB]) + (mu * patients[resistantB_drugB]) + (mu * patients[multiResistant_drugB]) - (3 * mu * patients[susceptible_drugB]);
        clearance =  (pharmacoGrowth[susceptible_drugB] - immuneClearance) * patients[susceptible_drugB];
        infection =  0.5 * (patients[susceptible_drugB] + patients[susceptible_drugA]) * patients[uninfected] * rateOfInfection;
        superInfection = ((patients[susceptible_drugA] + patients[susceptible_drugB]) * rateOfSuperInfection * rateOfInfection * cost) * (patients[resistantA_drugB] + patients[resistantB_drugB] + patients[multiResistant_drugB]);

        deltaS = mutation + clearance + infection + superInfection ;

        return deltaS;
    }

    private double deltaResistantA_drugA() {
        double deltaResistantA = 0;
        double mutation, clearance, infection, superInfection;

        mutation = (mu * patients[susceptible_drugA]) + (mu * patients[resistantB_drugA]) + (mu * patients[multiResistant_drugA]) - (3 * mu * patients[resistantA_drugA]);
        clearance = (pharmacoGrowth[resistantA_drugA] - immuneClearance) * patients[resistantA_drugA];
        infection =  rateOfInfection * patients[uninfected] * ((0.5 * (patients[resistantA_drugB] + patients[resistantA_drugA])) * (1.0-cost)) ;
        superInfection = patients[resistantA_drugA] * rateOfSuperInfection * rateOfInfection * cost * (patients[susceptible_drugA] + patients[susceptible_drugB]);

        deltaResistantA = mutation + clearance + infection - superInfection;

        return deltaResistantA;
    }

    private double deltaResistantA_drugB() {
        double deltaResistantA = 0;
        double mutation, clearance, infection, superInfection;

        mutation = (mu * patients[susceptible_drugB]) + (mu * patients[resistantB_drugB]) + (mu * patients[multiResistant_drugB]) - (3 * mu * patients[resistantA_drugB]);
        clearance = (pharmacoGrowth[resistantA_drugB] - immuneClearance) * patients[resistantA_drugB];
        infection =  (0.5 * (patients[resistantA_drugA] + patients[resistantA_drugB]) * (1.0-cost)) * patients[uninfected] * rateOfInfection;
        superInfection = patients[resistantA_drugB] * rateOfSuperInfection * rateOfInfection * cost * (patients[susceptible_drugA] + patients[susceptible_drugB]);

        deltaResistantA = mutation + clearance + infection - superInfection;

        return deltaResistantA;
    }

    private double deltaResistantB_drugA() {
        double deltaResistantB = 0;
        double mutation, clearance, infection, superInfection;


        mutation = (mu * patients[susceptible_drugA]) + (mu * patients[resistantA_drugA]) + (mu * patients[multiResistant_drugA]) - (3 * mu * patients[resistantB_drugA]);
        clearance = (pharmacoGrowth[resistantB_drugA] - immuneClearance) * patients[resistantB_drugA];
        infection =  (0.5 * (patients[resistantB_drugA] + patients[resistantB_drugB])) * (1.0-cost) * patients[uninfected] * rateOfInfection;
        superInfection = patients[resistantB_drugA]  * rateOfSuperInfection * rateOfInfection * cost * (patients[susceptible_drugA] + patients[susceptible_drugB]);

        deltaResistantB = mutation + clearance + infection - superInfection;

        return deltaResistantB;
    }

    private double deltaResistantB_drugB() {
        double deltaResistantB = 0;
        double mutation, clearance, infection, superInfection;


        mutation = (mu * patients[susceptible_drugB]) + (mu * patients[resistantA_drugB]) + (mu * patients[multiResistant_drugB]) - (3 * mu * patients[resistantB_drugB]);
        clearance = (pharmacoGrowth[resistantB_drugB] - immuneClearance) * patients[resistantB_drugB];
        infection =  (0.5 * (patients[resistantB_drugA] + patients[resistantB_drugB])) * (1.0-cost) * patients[uninfected] * rateOfInfection;
        superInfection = patients[resistantB_drugB] * rateOfSuperInfection * rateOfInfection * cost * (patients[susceptible_drugA] + patients[susceptible_drugB]);

        deltaResistantB = mutation + clearance + infection - superInfection;

        return deltaResistantB;
    }

    private double deltaMultiResistant_drugA() {
        double deltaMultiResistant = 0;
        double mutation, clearance, infection, superInfection;


        mutation = (mu * patients[susceptible_drugA]) + (mu * patients[resistantA_drugA]) + (mu * patients[resistantB_drugA]) - (3 * mu * patients[multiResistant_drugA]);
        clearance = (pharmacoGrowth[multiResistant_drugA] - immuneClearance) * patients[multiResistant_drugA];
        infection =  (0.5 * (patients[multiResistant_drugB] + patients[multiResistant_drugA])) * (1.0-cost) * patients[uninfected] * rateOfInfection;
        superInfection = patients[multiResistant_drugA] * rateOfSuperInfection * rateOfInfection * cost * (patients[susceptible_drugA] + patients[susceptible_drugB]);

        deltaMultiResistant = mutation + clearance + infection - superInfection;

        return deltaMultiResistant;
    }

    private double deltaMultiResistant_drugB() {
        double deltaMultiResistant = 0;
        double mutation, clearance, infection, superInfection;


        mutation = (mu * patients[susceptible_drugB]) + (mu * patients[resistantA_drugB]) + (mu * patients[resistantB_drugB]) - (3 * mu * patients[multiResistant_drugB]);
        clearance = (pharmacoGrowth[multiResistant_drugB] - immuneClearance) * patients[multiResistant_drugB];
        infection =  0.5 * (patients[multiResistant_drugA] + patients[multiResistant_drugB]) * (1.0-cost) * patients[uninfected] * rateOfInfection;
        superInfection = patients[multiResistant_drugB] * rateOfSuperInfection * rateOfInfection * cost * (patients[susceptible_drugA] + patients[susceptible_drugB]);

        deltaMultiResistant = mutation + clearance + infection - superInfection;

        return deltaMultiResistant;
    }

    private double deltaUninfected() {
        double deltaUninfected = 0;

        // recovery
        deltaUninfected -= (pharmacoGrowth[susceptible_drugA] - immuneClearance) * patients[susceptible_drugA];
        deltaUninfected -= (pharmacoGrowth[susceptible_drugB] - immuneClearance) * patients[susceptible_drugB];

        deltaUninfected -= (pharmacoGrowth[resistantA_drugA] - immuneClearance) * patients[resistantA_drugA];
        deltaUninfected -= (pharmacoGrowth[resistantA_drugB] - immuneClearance) * patients[resistantA_drugB];

        deltaUninfected -= (pharmacoGrowth[resistantB_drugA] - immuneClearance) * patients[resistantB_drugA];
        deltaUninfected -= (pharmacoGrowth[resistantB_drugB] - immuneClearance) * patients[resistantB_drugB];

        deltaUninfected -= (pharmacoGrowth[multiResistant_drugA] - immuneClearance) * patients[multiResistant_drugA];
        deltaUninfected -= (pharmacoGrowth[multiResistant_drugB] - immuneClearance) * patients[multiResistant_drugB];


        // infection
        deltaUninfected -= ((patients[susceptible_drugA] + patients[susceptible_drugB]) * rateOfInfection * patients[uninfected]);

        deltaUninfected -= ((patients[resistantA_drugA] + patients[resistantA_drugB]) * rateOfInfection * (1.0 - cost) * patients[uninfected]) ;

        deltaUninfected -= ((patients[resistantB_drugA] + patients[resistantB_drugB]) * rateOfInfection * (1.0 - cost) * patients[uninfected]) ;

        deltaUninfected -= ((patients[multiResistant_drugA] + patients[multiResistant_drugB]) * rateOfInfection * (1.0 - cost) * patients[uninfected]);

        return deltaUninfected;

    }

    private void findMeans() {
        double[] values = new double[9];
        int counter = 0;
        for (int i = burn; i < frequencies.size(); i++) {

            values[uninfected] += frequencies.get(i)[uninfected];

            values[susceptible_drugA] += frequencies.get(i)[susceptible_drugA];
            values[susceptible_drugB] += frequencies.get(i)[susceptible_drugB];

            values[resistantA_drugA] += frequencies.get(i)[resistantA_drugA];
            values[resistantA_drugB] += frequencies.get(i)[resistantA_drugB];

            values[resistantB_drugA] += frequencies.get(i)[resistantB_drugA];
            values[resistantB_drugB] += frequencies.get(i)[resistantB_drugB];

            values[multiResistant_drugA] += frequencies.get(i)[multiResistant_drugA];
            values[multiResistant_drugB] += frequencies.get(i)[multiResistant_drugB];

            counter++;

        }

        System.out.println((1.0 - (SimulationSettings.getInstance().getTradeoff()/2.0)) + "\t" + values[uninfected] / counter + "\t" + values[susceptible_drugA] / counter + "\t" + values[susceptible_drugB] / counter + "\t" + values[resistantA_drugA] / counter + "\t" + values[resistantA_drugB] / counter + "\t" + values[resistantB_drugA] / counter + "\t" + values[resistantB_drugB] / counter + "\t" + values[multiResistant_drugA] / counter + "\t" + values[multiResistant_drugB] / counter);
    }


}
