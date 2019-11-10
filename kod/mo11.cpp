#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

class DyskretyzacjaRozwiazanie
{
private:
    double xmin, xmax, tmin, tmax;
    double D, h, lambda, dt;
    int wiersz, kolumna;
    double PI = 3.14159265359;
public:
    DyskretyzacjaRozwiazanie(double lambda) :
        lambda(lambda), xmin(0), xmax(1.0),
        tmin(0), tmax(0.5), D(1.0), h(0.1)
    {
        dt = (lambda*h*h)/D;
        wiersz = static_cast<int>((tmax - tmin)/dt + 1);
        kolumna= static_cast<int>((xmax - xmin)/h + 1);
    }

    std::vector<std::vector<double>> Analityczna(std::vector<std::vector<double>> U)
    {
        for(int i = 1; i<U.size(); i++)
        {
             for(int j=1; j<U[0].size()-1; j++)
            {
                U[i][j] = (1 - exp(-(PI*PI)*this->D*dt*i))*sin(PI*j*h);
            }
        }
        return U;

    }

    std::vector<std::vector<double>> KlasycznaMetodaBezposrednia(std::vector<std::vector<double>> U){
        U=wypelnijMacierzWarunkiemBrzegowym(U);
        U=wypelnijMacierzWarunkiemPoczatkowym(U);

        for(int i = 1; i < U.size(); i++){
            for(int j = 1; j < U[0].size()-1; j++){
                 U[i][j]=lambda*U[i-1][j-1]+(1-2*lambda)*U[i-1][j]+lambda*U[i-1][j+1]+dt*(PI*PI)*sin(PI*j*h);
            }
        }
        return U;
    }


    int wyborCzesciowy(std::vector<std::vector<double>>& macierz, int j, int *indeksy){
        int w;
        for(int i = j + 1; i < macierz.size(); i++) {
            if(fabs(macierz.at(indeksy[i]).at(j)) < fabs(macierz.at(indeksy[i + 1]).at(j))) {
                w = indeksy[i + 1];
            }
            else{
                w = indeksy[i];
            }
        }
        return w;
    }

    std::vector<std::vector<double>> dekompozycjaLU(std::vector<std::vector<double>> macierz){

        int rozmiar = macierz.size();
        int* indeksy = new int[rozmiar]; //wektor przechowywuj¹cy kolejnoœæ wierszy
        double element;

        for(int i = 0; i < rozmiar; i++)
            indeksy[i] = i; //numeracja wierszy macierzy po kolei

        //k oznacza k etap
        for(int k = 1; k < rozmiar; k++) {
            element = macierz.at(indeksy[k - 1]).at(k - 1);

            if(element == 0.0) {
                //bierzemy biezacy (dla nas to bedzie 1 czyli kolejny) wiersz i iterujemy do max rozmiaru
                int w = wyborCzesciowy(macierz, indeksy[k], indeksy);
                indeksy[w] = indeksy[k];
                indeksy[k] = w;
                element = macierz[indeksy[k - 1]][k - 1];
            }
            //rozpoczecie iterowania dla tablicy zaczynajacej sie od przekatnej
            for(int i = k; i < rozmiar; i++) {
                //zaczynamy od 'wiersza nizej' bo i=k
                double mnoznik = macierz[indeksy[i]][k - 1] / element;
                macierz[indeksy[i]][k - 1]= mnoznik;
                for(int j = k; j < rozmiar; j++) {
                    macierz[indeksy[i]][j] -= macierz[indeksy[k - 1]][j] * mnoznik;
                }
            }
        }
        return macierz;
    }

    std::vector<double> rozwiaz(std::vector<std::vector<double>> macierz, std::vector<double> wektor){
        int rozmiar = macierz.size();
        std::vector<double> x(rozmiar, 0);
        std::vector<double> y(rozmiar, 0);
        y[0] = wektor[0]; //dla 0 0 ze wzoru
        for(int i = 1; i < rozmiar; i++) {
            y[i] = wektor[i];
            for(int j = 0; j < i; j++) {
                y[i] = y[i] - macierz[i][j] * y[j];
            }
        }
        x[rozmiar - 1] = y[rozmiar - 1] / macierz[rozmiar - 1][rozmiar - 1];
        for(int i = rozmiar - 1; i >= 0; i--) {
            x[i] = y[i];
            for(int j = i + 1; j < rozmiar; j++) {
                //od przekatnej do konca
                x[i] = x[i] - macierz[i][j] * x[j];
            }
            //koncowe podzielenie przez wartosc na przekatnej
            x[i] = x[i] / macierz[i][i];


        }
        return x;
    }


    std::vector<std::vector<double>> CrankNicolson(std::vector<std::vector<double>> U){
        U=wypelnijMacierzWarunkiemBrzegowym(U);
        U=wypelnijMacierzWarunkiemPoczatkowym(U);
        std::vector<std::vector<double>> A(kolumna, std::vector<double>(kolumna, 0));
        std::vector<double> B;

        A.at(0).at(0) = 1;

        for(int i = 1; i<kolumna-1; i++){
            A.at(i).at(i-1) = lambda/2.0;
            A.at(i).at(i) = -(1+lambda);
            A.at(i).at(i+1) = lambda/2.0;
        }
        A.at(kolumna-1).at(kolumna-1) = 1;
        std::vector<std::vector<double>> LU = dekompozycjaLU(A);
        for(int i = 1; i<wiersz; i++){
            B.clear();
            B.push_back(0);
            for(int j = 1; j<kolumna-1; j++){
                B.push_back(-((lambda/2.0)*U.at(i-1).at(j-1)+(1-lambda)*U.at(i-1).at(j)+(lambda/2.0)*U.at(i-1).at(j+1))-dt*(PI*PI)*sin(PI*j*h));
            }
            B.push_back(0);
            U.at(i) = rozwiaz(LU, B);
        }
        return U;
    }



    std::vector<std::vector<double>> liczMacierzBledu(std::vector<std::vector<double>> U, std::vector<std::vector<double>> Uanali)
    {
        std::vector<std::vector<double>> blad(U.size(), std::vector<double>(U[0].size(), 0));
        for(int i = 0; i<U.size(); i++){
            for(int j = 0; j<U[0].size();j++){
                blad.at(i).at(j) = fabs(U[i][j] - Uanali[i][j]);
            }
        }
        return blad;
    }

    std::vector<double> liczWektorBledowMaksOdT(std::vector<std::vector<double>> macierzBledow)
    {
        std::vector<double> blad(macierzBledow.size(), 0);
        double max;
        for(int i = 0; i<macierzBledow.size(); i++){
            blad[i] = fabs(macierzBledow[i][0]);
            for(int j = 0; j<macierzBledow[0].size();j++){
                if(fabs(blad[i])<fabs(macierzBledow[i][j]))
                    blad[i]=fabs(macierzBledow[i][j]);
            }
        }
        return blad;
    }

    std::vector<std::vector<double>> wypelnijMacierzWarunkiemBrzegowym(std::vector<std::vector<double>> macierz){
        for(int i = 0; i<macierz.size();i++){
            macierz[i][0] = 0;
            macierz[i][macierz[0].size()-1] = 0;}
        return macierz;}


    std::vector<std::vector<double>> wypelnijMacierzWarunkiemPoczatkowym(std::vector<std::vector<double>> macierz){
        for(int i = 0; i<macierz[0].size();i++){
            macierz[0][i] = 0;}
        return macierz;}


    void zapiszRozwiazanieTransponowaneDoPliku(std::fstream& plik, std::vector<std::vector<double >> U){
        std::vector<std::vector<double >> UT = transponujMacierz(U);
        for(int i=0;i<UT.size();i++){
            plik << this->h*i << " ";
            for(int j = 0;j<UT[0].size();j++){
                plik << UT[i][j] << " ";}
        plik << std::endl;}}


    void zapiszRozwiazanieDoPliku(std::fstream& plik, std::vector<std::vector<double >> U){
        for(int i=0;i<U.size();i++){
            for(int j = 0;j<U[0].size();j++){
                plik << U[i][j] << " ";}
        plik << std::endl;}}

    void zapiszMaxBladDoPliku(std::fstream& plik, std::vector<double> wektor){
        for(int i = 0; i<wektor.size();i++){
            plik << dt*i << " " << wektor.at(i) << " " <<  log10(dt*i) << " " << log10(wektor.at(i)) << std::endl;}}

    void zapiszMaxBladTMaksDoPliku(std::fstream& plik, std::vector<double> wektor){
        double k = 0.01;
        for(int i = 0; i<wektor.size();i++){
            plik << k << " " << wektor.at(i) << " " << log10(k) << " " << log10(wektor.at(i)) << std::endl;
            k+=0.01;}}

    std::vector<std::vector<double >> transponujMacierz(std::vector<std::vector<double >> U){
        std::vector<std::vector<double>> UT(U[0].size(), std::vector<double>(U.size(), 0));
        for(int i = 0; i < U.size(); i++)
            for(int j = 0; j < U[0].size(); j++) UT[j][i] = U[i][j];
        return UT;}


    void rysujMacierz(std::vector<std::vector<double>> macierz){
        for(int i = 0; i<macierz.size();i++){
            for(int j = 0; j<macierz[0].size();j++){
                std::cout << macierz.at(i).at(j);}
            std::cout << std::endl;}}


    void rysujWektor(std::vector<double> wektor){
        for(int i = 0; i<wektor.size();i++){
            std::cout << wektor.at(i);}}

    std::vector<std::vector<double>> dajMacierz(){
        std::vector<std::vector<double>> U(wiersz, std::vector<double>(kolumna, 0));
        return U;}



};


int main()
{
    DyskretyzacjaRozwiazanie KMB(0.4);
    std::vector<std::vector<double>> U1 = KMB.dajMacierz();
    std::fstream fakmb("analityczneKMB.txt", std::ios::out);
    std::fstream fbkmb("KMB.txt", std::ios::out);
    std::fstream fbladkmb("macierzbleduKMB.txt", std::ios::out);
    std::fstream fmaxbladkmb("wektormaxbleduKMB.txt", std::ios::out);

    std::vector<std::vector<double>> AKMB;
    std::vector<std::vector<double>> BKMB;
    std::vector<std::vector<double>> bladKMB;
    std::vector<double> maksBladKMB;

    AKMB = KMB.Analityczna(U1);
    BKMB = KMB.KlasycznaMetodaBezposrednia(U1);
    bladKMB = KMB.liczMacierzBledu(BKMB,AKMB);
    maksBladKMB = KMB.liczWektorBledowMaksOdT(bladKMB);

    KMB.zapiszRozwiazanieTransponowaneDoPliku(fakmb,AKMB);
    KMB.zapiszRozwiazanieTransponowaneDoPliku(fbkmb,BKMB);
    KMB.zapiszRozwiazanieTransponowaneDoPliku(fbladkmb,bladKMB);
    KMB.zapiszMaxBladDoPliku(fmaxbladkmb, maksBladKMB);

    //

    DyskretyzacjaRozwiazanie CN(1.0);
    std::vector<std::vector<double>> U2 = CN.dajMacierz();

    std::fstream facn("analityczneCN.txt", std::ios::out);
    std::fstream fbcn("CN.txt", std::ios::out);
    std::fstream fbladcn("macierzbleduCN.txt", std::ios::out);
    std::fstream fmaxbladcn("wektormaxbleduCN.txt", std::ios::out);

    std::vector<std::vector<double>> ACN;
    std::vector<std::vector<double>> BCN;
    std::vector<std::vector<double>> bladCN;
    std::vector<double> maksBladCN;

    ACN = CN.Analityczna(U2);
    BCN = CN.CrankNicolson(U2);
    bladCN = CN.liczMacierzBledu(BCN,ACN);
    maksBladCN = CN.liczWektorBledowMaksOdT(bladCN);

    CN.zapiszRozwiazanieTransponowaneDoPliku(facn,ACN);
    CN.zapiszRozwiazanieTransponowaneDoPliku(fbcn,BCN);
    CN.zapiszRozwiazanieTransponowaneDoPliku(fbladcn,bladCN);
    CN.zapiszMaxBladDoPliku(fmaxbladcn, maksBladCN);

}
