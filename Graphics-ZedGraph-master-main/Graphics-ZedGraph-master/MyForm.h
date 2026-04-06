#pragma once
#include <math.h>
#include <iostream>
#include <vector>

/* --- Test Task --- */
struct DataTestTask {
	double x;
	double V_full;
	double V_half2;
	double Vfull_Vhalf2;
	double OLP;
	double h;
	int C1;
	int C2;
	double U;
	double U_V;
};
double variant = 7;
double f_test_task(double x, double u) { 
	double l = ((int)variant % 2 ? -1 : 1) * (variant / 2.0);
	return l * u;
}
static double f_test_pervoobr(double x) {
	return exp(-7.0 / 2.0 * x);
}
double RK4_step_test_task(double V, double x, double h) {
	double k1 = f_test_task(x, V);
	double k2 = f_test_task(x + h / 2.0, V + h * k1 / 2.0);
	double k3 = f_test_task(x + h / 2.0, V + h * k2 / 2.0);
	double k4 = f_test_task(x + h, V + h * k3);

	return V + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
}
std::vector<DataTestTask> RK4_method_fixed_step_test_task(double x0, double U0, double h, double b, int Nmax) {
	std::vector<DataTestTask> Data;

	double x_n = x0;
	double V_n = U0;
	int iters = 0;

	DataTestTask CurrStepData;

	CurrStepData.x = x_n;
	CurrStepData.V_full = V_n;
	CurrStepData.V_half2 = V_n;
	CurrStepData.Vfull_Vhalf2 = V_n - V_n;
	CurrStepData.OLP = 16 * (V_n - V_n) / 15.0;
	CurrStepData.h = h;
	CurrStepData.C1 = 0;
	CurrStepData.C2 = 0;
	CurrStepData.U = f_test_pervoobr(x_n);
	CurrStepData.U_V = fabs(f_test_pervoobr(x_n) - V_n);

	Data.push_back(CurrStepData);


	while ((h > 0 && x_n < b) || (h < 0 && x_n > b)){
		DataTestTask CurrStepData;
		iters++;

		if (iters > Nmax) {
			std::cout << "Too many iterations\n";
			return Data;
		}

		if ((h > 0 && x_n + h > b) || (h < 0 && x_n + h < b)) {
			h = b - x_n;
		}

		double V_full = RK4_step_test_task(V_n, x_n, h);
		

		// Подсчет двойного шага для таблицы
		double V_half = RK4_step_test_task(V_n, x_n, h / 2.0);
		double V_half2 = RK4_step_test_task(V_half, x_n + h / 2.0, h / 2.0);

		V_n = V_full;
		x_n += h;

		CurrStepData.x = x_n;
		CurrStepData.V_full = V_n;
		CurrStepData.V_half2 = V_half2;
		CurrStepData.Vfull_Vhalf2 = V_n - V_half2;
		CurrStepData.OLP = 16 * (V_n - V_half2) / 15.0;
		CurrStepData.h = h;
		CurrStepData.C1 = 0;
		CurrStepData.C2 = 0;
		CurrStepData.U = f_test_pervoobr(x_n);
		CurrStepData.U_V = fabs(f_test_pervoobr(x_n) - V_n);

		Data.push_back(CurrStepData);

		if (fabs(x_n - b) < 1e-12) {
			break;
		}
	}

	return Data;
}
std::vector<DataTestTask> RK4_method_addaptive_step_test_task(double x0, double U0, double h0, double b, int Nmax, double Eps){
	std::vector<DataTestTask> Data;

	double x_n = x0;
	double V_n = U0;
	double h_n = h0;

	int iters = 0;

	double C1 = 0;
	double C2 = 0;

	DataTestTask CurrStepData;

	CurrStepData.x = x_n;
	CurrStepData.V_full = V_n;
	CurrStepData.V_half2 = V_n;
	CurrStepData.Vfull_Vhalf2 = V_n - V_n;
	CurrStepData.OLP = 16 * (V_n - V_n) / 15.0;
	CurrStepData.h = h_n;
	CurrStepData.C1 = C1;
	CurrStepData.C2 = C2;
	CurrStepData.U = f_test_pervoobr(x_n);
	CurrStepData.U_V = fabs(f_test_pervoobr(x_n) - V_n);

	Data.push_back(CurrStepData);

	while ((h_n > 0 && x_n < b) || (h_n < 0 && x_n > b)) {
		DataTestTask CurrStepData;

		iters++;

		if (iters > Nmax) {
			std::cout << "Too many iterations\n";
			return Data;
		}

		if ((h_n > 0 && x_n + h_n > b) || (h_n < 0 && x_n + h_n < b)) {
			h_n = b - x_n;
		}

		// полный шаг
		double V_full = RK4_step_test_task(V_n, x_n, h_n);

		// два полушага
		double V_half = RK4_step_test_task(V_n, x_n, h_n / 2.0);
		double V_half2 = RK4_step_test_task(V_half, x_n + h_n / 2.0, h_n / 2.0);

		// оценка локальной погрешности
		double S = (V_half2 - V_full) / 15.0;

		double h_temp = h_n;

		if (fabs(S) > Eps) {
			h_n *= 0.5;
			C1++;
			continue;
		}

		if (fabs(S) < Eps / 32.0) {
			h_n *= 2.0;
			C2++;
		}

		x_n += h_temp;
		V_n = V_half2;

		CurrStepData.x = x_n;
		CurrStepData.V_full = V_full;
		CurrStepData.V_half2 = V_half2;
		CurrStepData.Vfull_Vhalf2 = V_full - V_half2;
		CurrStepData.OLP = 16.0 * (V_full - V_half2) / 15.0;
		CurrStepData.h = h_n;
		CurrStepData.C1 = C1;
		CurrStepData.C2 = C2;
		CurrStepData.U = f_test_pervoobr(x_n);
		CurrStepData.U_V = fabs(f_test_pervoobr(x_n) - V_n);

		Data.push_back(CurrStepData);

		if (fabs(x_n - b) < 1e-12) {
			break;
		}
	}

	return Data;
}
/* --- Main Task --- */
struct DataMainTask {
	double x;
	double V_full;
	double V_half2;
	double Vfull_Vhalf2;
	double OLP;
	double h;
	int C1;
	int C2;
	double proizv;
};
double f1_main_task(double x, double U1, double U2) {
	return U2;
}
double f2_main_task(double x, double U1, double U2) {
	double k = 2.0;
	double k_ = 2.0;
	double c = 0.15;
	double m = 0.01; 

	return -(c / m) * U2 - (k / m) * U1 - (k_ / m) * pow(U1, 3);
}
struct State
{
	double U1;
	double U2;
};

State RK4_step_main_task(State S, double x, double h) {
	State k1;
	k1.U1 = f1_main_task(x, S.U1, S.U2);
	k1.U2 = f2_main_task(x, S.U1, S.U2);
	State k2;
	k2.U1 = f1_main_task(x + h / 2, S.U1 + h * k1.U1 / 2, S.U2 + h * k1.U2 / 2);
	k2.U2 = f2_main_task(x + h / 2, S.U1 + h * k1.U1 / 2, S.U2 + h * k1.U2 / 2);
	State k3;
	k3.U1 = f1_main_task(x + h / 2, S.U1 + h * k2.U1 / 2, S.U2 + h * k2.U2 / 2);
	k3.U2 = f2_main_task(x + h / 2, S.U1 + h * k2.U1 / 2, S.U2 + h * k2.U2 / 2);
	State k4;
	k4.U1 = f1_main_task(x + h, S.U1 + h * k3.U1, S.U2 + h * k3.U2);
	k4.U2 = f2_main_task(x + h, S.U1 + h * k3.U1, S.U2 + h * k3.U2);
	State result;
	result.U1 = S.U1 + h * (k1.U1 + 2 * k2.U1 + 2 * k3.U1 + k4.U1) / 6.0;
	result.U2 = S.U2 + h * (k1.U2 + 2 * k2.U2 + 2 * k3.U2 + k4.U2) / 6.0;

	return result;
}
std::vector<DataMainTask> RK4_method_fixed_step_main_task(double x0, double U1_0, double U2_0, double h, double b, int Nmax) {
	std::vector<DataMainTask> Data;

	double x_n = x0;
	double V1_n = U1_0;
	double V2_n = U2_0;
	int iters = 0;

	State V_n;
	V_n.U1 = V1_n;
	V_n.U2 = V2_n;

	DataMainTask CurrStepData;

	CurrStepData.x = x_n;
	CurrStepData.V_full = V_n.U1; // ?

	CurrStepData.V_half2 = V_n.U1; // ?
	CurrStepData.Vfull_Vhalf2 = V_n.U1 - V_n.U1; // ?
	CurrStepData.OLP = (V_n.U1 - V_n.U1) / 15.0; // ?
	CurrStepData.h = h;
	CurrStepData.C1 = 0;
	CurrStepData.C2 = 0;
	CurrStepData.proizv = U2_0;
	Data.push_back(CurrStepData);


	while ((h > 0 && x_n < b) || (h < 0 && x_n > b)) {
		DataMainTask CurrStepData;
		iters++;

		if (iters > Nmax) {
			std::cout << "Too many iterations\n";
			return Data;
		}

		if ((h > 0 && x_n + h > b) || (h < 0 && x_n + h < b)) {
			h = b - x_n;
		}

		State V_n1;
		V_n1 = RK4_step_main_task(V_n, x_n, h);

		// два полушага для таблицы
		State V_half;
		V_half = RK4_step_main_task(V_n, x_n, h / 2.0);
		State V_half2;
		V_half2 = RK4_step_main_task(V_half, x_n + h / 2.0, h / 2.0);

		CurrStepData.x = x_n;
		CurrStepData.V_full = V_n1.U1; // ?
		CurrStepData.V_half2 = V_half2.U1; // ?
		CurrStepData.Vfull_Vhalf2 = V_n1.U1 - V_half2.U1; // ?
		CurrStepData.OLP = (V_n1.U1 - V_half2.U1) / 15.0; // ?
		CurrStepData.h = h;
		CurrStepData.C1 = 0;
		CurrStepData.C2 = 0;
		CurrStepData.proizv = V_n1.U2;


		V_n = V_n1;
		x_n += h;

		Data.push_back(CurrStepData);

		if (fabs(x_n - b) < 1e-12) {
			break;
		}
	}

	return Data;
}
std::vector<DataMainTask> RK4_method_addaptive_step_main_task(double x0, double U1_0, double U2_0, double h, double b, int Nmax, double Eps) {
	std::vector<DataMainTask> Data;

	double x_n = x0;
	double V1_n = U1_0;
	double V2_n = U2_0;
	double h_n = h;
	int iters = 0;

	double C1 = 0;
	double C2 = 0;

	State V_n;
	V_n.U1 = V1_n;
	V_n.U2 = V2_n;

	DataMainTask CurrStepData;

	CurrStepData.x = x_n;
	CurrStepData.V_full = V_n.U1; // ?
	CurrStepData.V_half2 = V_n.U1; // ?
	CurrStepData.Vfull_Vhalf2 = V_n.U1 - V_n.U1; // ?
	CurrStepData.OLP = (V_n.U1 - V_n.U1) / 15.0; // ?
	CurrStepData.h = h;
	CurrStepData.C1 = 0;
	CurrStepData.C2 = 0;
	CurrStepData.proizv = U2_0;


	Data.push_back(CurrStepData);

	while ((h_n > 0 && x_n < b) || (h_n < 0 && x_n > b)) {
		DataMainTask CurrStepData;
		iters++;

		if (iters > Nmax) {
			std::cout << "Too many iterations\n";
			return Data;
		}

		if ((h_n > 0 && x_n + h_n > b) || (h_n < 0 && x_n + h_n < b)) {
			h_n = b - x_n;
		}

		// полный шаг
		State V_full;
		V_full = RK4_step_main_task(V_n, x_n, h_n);

		// два полушага
		State V_half;
		V_half = RK4_step_main_task(V_n, x_n, h_n / 2.0);
		State V_half2;
		V_half2 = RK4_step_main_task(V_half, x_n + h_n / 2.0, h_n / 2.0);

		// оценка локальной погрешности
		double S1 = (V_half2.U1 - V_full.U1) / 15.0;
		double S2 = (V_half2.U2 - V_full.U2) / 15.0;

		double h_temp = h_n;

		if (fabs(S1) > Eps || fabs(S2) > Eps) {
			h_n *= 0.5;
			C1++;
			continue;
		}

		if ((fabs(S1) < Eps / 32.0) && (fabs(S2) < Eps / 32.0)) {
			h_n *= 2.0;
			C2++;
		}

		V_n = V_half2;
		x_n += h_temp;

		CurrStepData.x = x_n;
		CurrStepData.V_full = V_full.U1; // ?
		CurrStepData.V_half2 = V_half2.U1; // ?
		CurrStepData.Vfull_Vhalf2 = V_full.U1 - V_half2.U1; // ?
		CurrStepData.OLP = (V_full.U1 - V_half2.U1) / 15.0; // ?
		CurrStepData.h = h_n;
		CurrStepData.C1 = C1;
		CurrStepData.C2 = C2;
		CurrStepData.proizv = V_full.U2;


		Data.push_back(CurrStepData);

		if (fabs(x_n - b) < 1e-12) {
			break;
		}
	}

	return Data;
}

namespace Graph {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace ZedGraph;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: ZedGraph::ZedGraphControl^  zedGraphControl1;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::DataGridView^  dataGridView1;



	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::TextBox^  textBox3;
	private: System::Windows::Forms::Button^  button2;
	private: System::Windows::Forms::TextBox^  textBox4;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TextBox^  textBox5;
	private: System::Windows::Forms::Label^  label5;










	private: System::Windows::Forms::TextBox^ textBox6;
	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::TabControl^ tabControl1;
	private: System::Windows::Forms::TabPage^ TestTask;
	private: System::Windows::Forms::TabPage^ MainTask;
	private: System::Windows::Forms::CheckBox^ checkBox1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ Iter;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ X;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ V_full;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ V_half2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ V_full_V_half2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ OLP;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ h;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ C1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ C2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ U;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ U_V;
	private: ZedGraph::ZedGraphControl^ zedGraphControl2;
	private: System::Windows::Forms::DataGridView^ dataGridView2;











	private: System::Windows::Forms::CheckBox^ checkBox2;
	private: System::Windows::Forms::Button^ button3;
	private: System::Windows::Forms::Label^ label7;
	private: System::Windows::Forms::TextBox^ textBox7;
	private: System::Windows::Forms::Button^ button4;
	private: System::Windows::Forms::TextBox^ textBox8;
	private: System::Windows::Forms::Label^ label8;
	private: System::Windows::Forms::TextBox^ textBox9;
	private: System::Windows::Forms::Label^ label9;
	private: System::Windows::Forms::Label^ label10;
	private: System::Windows::Forms::TextBox^ textBox10;
	private: System::Windows::Forms::TextBox^ textBox11;
	private: System::Windows::Forms::Label^ label11;
	private: System::Windows::Forms::Label^ label12;
	private: System::Windows::Forms::TextBox^ textBox12;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn2;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn3;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn4;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn5;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn6;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn7;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn8;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ dataGridViewTextBoxColumn9;
private: System::Windows::Forms::TextBox^ textBox13;
private: System::Windows::Forms::Label^ label13;
private: System::Windows::Forms::TabControl^ tabControl2;
private: System::Windows::Forms::TabPage^ ux;
private: System::Windows::Forms::TabPage^ uproizvx;
private: ZedGraph::ZedGraphControl^ zedGraphControl3;
private: System::Windows::Forms::TabPage^ phase;
private: ZedGraph::ZedGraphControl^ zedGraphControl4;




























	protected:
	private: System::ComponentModel::IContainer^  components;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->zedGraphControl1 = (gcnew ZedGraph::ZedGraphControl());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->Iter = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->X = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->V_full = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->V_half2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->V_full_V_half2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->OLP = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->h = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->C1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->C2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->U = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->U_V = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->textBox5 = (gcnew System::Windows::Forms::TextBox());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->textBox6 = (gcnew System::Windows::Forms::TextBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->tabControl1 = (gcnew System::Windows::Forms::TabControl());
			this->TestTask = (gcnew System::Windows::Forms::TabPage());
			this->checkBox1 = (gcnew System::Windows::Forms::CheckBox());
			this->MainTask = (gcnew System::Windows::Forms::TabPage());
			this->textBox13 = (gcnew System::Windows::Forms::TextBox());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->dataGridView2 = (gcnew System::Windows::Forms::DataGridView());
			this->dataGridViewTextBoxColumn1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn3 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn4 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn5 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn6 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn7 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn8 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->dataGridViewTextBoxColumn9 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->checkBox2 = (gcnew System::Windows::Forms::CheckBox());
			this->button3 = (gcnew System::Windows::Forms::Button());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->textBox7 = (gcnew System::Windows::Forms::TextBox());
			this->button4 = (gcnew System::Windows::Forms::Button());
			this->textBox8 = (gcnew System::Windows::Forms::TextBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->textBox9 = (gcnew System::Windows::Forms::TextBox());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->textBox10 = (gcnew System::Windows::Forms::TextBox());
			this->textBox11 = (gcnew System::Windows::Forms::TextBox());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->textBox12 = (gcnew System::Windows::Forms::TextBox());
			this->zedGraphControl2 = (gcnew ZedGraph::ZedGraphControl());
			this->tabControl2 = (gcnew System::Windows::Forms::TabControl());
			this->ux = (gcnew System::Windows::Forms::TabPage());
			this->uproizvx = (gcnew System::Windows::Forms::TabPage());
			this->phase = (gcnew System::Windows::Forms::TabPage());
			this->zedGraphControl3 = (gcnew ZedGraph::ZedGraphControl());
			this->zedGraphControl4 = (gcnew ZedGraph::ZedGraphControl());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->tabControl1->SuspendLayout();
			this->TestTask->SuspendLayout();
			this->MainTask->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView2))->BeginInit();
			this->tabControl2->SuspendLayout();
			this->ux->SuspendLayout();
			this->uproizvx->SuspendLayout();
			this->phase->SuspendLayout();
			this->SuspendLayout();
			// 
			// zedGraphControl1
			// 
			this->zedGraphControl1->Location = System::Drawing::Point(0, 0);
			this->zedGraphControl1->Name = L"zedGraphControl1";
			this->zedGraphControl1->ScrollGrace = 0;
			this->zedGraphControl1->ScrollMaxX = 0;
			this->zedGraphControl1->ScrollMaxY = 0;
			this->zedGraphControl1->ScrollMaxY2 = 0;
			this->zedGraphControl1->ScrollMinX = 0;
			this->zedGraphControl1->ScrollMinY = 0;
			this->zedGraphControl1->ScrollMinY2 = 0;
			this->zedGraphControl1->Size = System::Drawing::Size(613, 375);
			this->zedGraphControl1->TabIndex = 0;
			this->zedGraphControl1->Load += gcnew System::EventHandler(this, &MyForm::zedGraphControl1_Load);
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(538, 393);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(70, 29);
			this->button1->TabIndex = 1;
			this->button1->Text = L"Draw";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// dataGridView1
			// 
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(11) {
				this->Iter, this->X,
					this->V_full, this->V_half2, this->V_full_V_half2, this->OLP, this->h, this->C1, this->C2, this->U, this->U_V
			});
			this->dataGridView1->Location = System::Drawing::Point(614, 0);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->RowHeadersVisible = false;
			this->dataGridView1->Size = System::Drawing::Size(957, 375);
			this->dataGridView1->TabIndex = 2;
			this->dataGridView1->CellContentClick += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MyForm::dataGridView1_CellContentClick);
			// 
			// Iter
			// 
			this->Iter->HeaderText = L"Iter";
			this->Iter->Name = L"Iter";
			this->Iter->Width = 40;
			// 
			// X
			// 
			this->X->HeaderText = L"X";
			this->X->Name = L"X";
			this->X->ReadOnly = true;
			this->X->Width = 70;
			// 
			// V_full
			// 
			this->V_full->HeaderText = L"V_full";
			this->V_full->Name = L"V_full";
			this->V_full->ReadOnly = true;
			// 
			// V_half2
			// 
			this->V_half2->HeaderText = L"V_half2";
			this->V_half2->Name = L"V_half2";
			this->V_half2->ReadOnly = true;
			// 
			// V_full_V_half2
			// 
			this->V_full_V_half2->HeaderText = L"V_full - V_half2";
			this->V_full_V_half2->Name = L"V_full_V_half2";
			// 
			// OLP
			// 
			this->OLP->HeaderText = L"ОЛП";
			this->OLP->Name = L"OLP";
			// 
			// h
			// 
			this->h->HeaderText = L"h";
			this->h->Name = L"h";
			this->h->Width = 70;
			// 
			// C1
			// 
			this->C1->HeaderText = L"C1";
			this->C1->Name = L"C1";
			this->C1->Width = 40;
			// 
			// C2
			// 
			this->C2->HeaderText = L"C2";
			this->C2->Name = L"C2";
			this->C2->Width = 40;
			// 
			// U
			// 
			this->U->HeaderText = L"U";
			this->U->Name = L"U";
			// 
			// U_V
			// 
			this->U_V->HeaderText = L"|U - V|";
			this->U_V->Name = L"U_V";
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(6, 396);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(13, 13);
			this->label1->TabIndex = 3;
			this->label1->Text = L"a";
			this->label1->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(25, 393);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(48, 20);
			this->textBox1->TabIndex = 4;
			this->textBox1->Text = L"0";
			this->textBox1->TextChanged += gcnew System::EventHandler(this, &MyForm::textBox1_TextChanged);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(91, 396);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(13, 13);
			this->label2->TabIndex = 5;
			this->label2->Text = L"b";
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(110, 393);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(49, 20);
			this->textBox2->TabIndex = 6;
			this->textBox2->Text = L"1";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(256, 396);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(13, 13);
			this->label3->TabIndex = 7;
			this->label3->Text = L"h";
			// 
			// textBox3
			// 
			this->textBox3->Location = System::Drawing::Point(275, 393);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(61, 20);
			this->textBox3->TabIndex = 8;
			this->textBox3->Text = L"0,1";
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(614, 393);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(80, 29);
			this->button2->TabIndex = 9;
			this->button2->Text = L"Zoom";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// textBox4
			// 
			this->textBox4->Location = System::Drawing::Point(376, 394);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(49, 20);
			this->textBox4->TabIndex = 13;
			this->textBox4->Text = L"0,001";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(346, 397);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(25, 13);
			this->label4->TabIndex = 12;
			this->label4->Text = L"Eps";
			this->label4->Click += gcnew System::EventHandler(this, &MyForm::label4_Click);
			// 
			// textBox5
			// 
			this->textBox5->Location = System::Drawing::Point(195, 393);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(48, 20);
			this->textBox5->TabIndex = 11;
			this->textBox5->Text = L"1,0";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(168, 397);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(21, 13);
			this->label5->TabIndex = 10;
			this->label5->Text = L"U0";
			this->label5->Click += gcnew System::EventHandler(this, &MyForm::label5_Click);
			// 
			// textBox6
			// 
			this->textBox6->Location = System::Drawing::Point(462, 394);
			this->textBox6->Name = L"textBox6";
			this->textBox6->Size = System::Drawing::Size(49, 20);
			this->textBox6->TabIndex = 14;
			this->textBox6->Text = L"20";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Location = System::Drawing::Point(431, 398);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(34, 13);
			this->label6->TabIndex = 15;
			this->label6->Text = L"Nmax";
			this->label6->Click += gcnew System::EventHandler(this, &MyForm::label6_Click);
			// 
			// tabControl1
			// 
			this->tabControl1->Controls->Add(this->TestTask);
			this->tabControl1->Controls->Add(this->MainTask);
			this->tabControl1->Location = System::Drawing::Point(0, 0);
			this->tabControl1->Name = L"tabControl1";
			this->tabControl1->SelectedIndex = 0;
			this->tabControl1->Size = System::Drawing::Size(1594, 679);
			this->tabControl1->TabIndex = 16;
			// 
			// TestTask
			// 
			this->TestTask->Controls->Add(this->checkBox1);
			this->TestTask->Controls->Add(this->zedGraphControl1);
			this->TestTask->Controls->Add(this->dataGridView1);
			this->TestTask->Controls->Add(this->button2);
			this->TestTask->Controls->Add(this->label6);
			this->TestTask->Controls->Add(this->textBox1);
			this->TestTask->Controls->Add(this->button1);
			this->TestTask->Controls->Add(this->textBox6);
			this->TestTask->Controls->Add(this->label1);
			this->TestTask->Controls->Add(this->textBox4);
			this->TestTask->Controls->Add(this->label2);
			this->TestTask->Controls->Add(this->label4);
			this->TestTask->Controls->Add(this->textBox2);
			this->TestTask->Controls->Add(this->textBox5);
			this->TestTask->Controls->Add(this->label3);
			this->TestTask->Controls->Add(this->label5);
			this->TestTask->Controls->Add(this->textBox3);
			this->TestTask->Location = System::Drawing::Point(4, 22);
			this->TestTask->Name = L"TestTask";
			this->TestTask->Padding = System::Windows::Forms::Padding(3);
			this->TestTask->Size = System::Drawing::Size(1586, 653);
			this->TestTask->TabIndex = 0;
			this->TestTask->Text = L"Test Task";
			this->TestTask->UseVisualStyleBackColor = true;
			this->TestTask->Click += gcnew System::EventHandler(this, &MyForm::TestTask_Click);
			// 
			// checkBox1
			// 
			this->checkBox1->AutoSize = true;
			this->checkBox1->Location = System::Drawing::Point(256, 419);
			this->checkBox1->Name = L"checkBox1";
			this->checkBox1->Size = System::Drawing::Size(76, 17);
			this->checkBox1->TabIndex = 16;
			this->checkBox1->Text = L"Fixed Step";
			this->checkBox1->UseVisualStyleBackColor = true;
			this->checkBox1->CheckedChanged += gcnew System::EventHandler(this, &MyForm::checkBox1_CheckedChanged);
			// 
			// MainTask
			// 
			this->MainTask->Controls->Add(this->tabControl2);
			this->MainTask->Controls->Add(this->textBox13);
			this->MainTask->Controls->Add(this->label13);
			this->MainTask->Controls->Add(this->dataGridView2);
			this->MainTask->Controls->Add(this->checkBox2);
			this->MainTask->Controls->Add(this->button3);
			this->MainTask->Controls->Add(this->label7);
			this->MainTask->Controls->Add(this->textBox7);
			this->MainTask->Controls->Add(this->button4);
			this->MainTask->Controls->Add(this->textBox8);
			this->MainTask->Controls->Add(this->label8);
			this->MainTask->Controls->Add(this->textBox9);
			this->MainTask->Controls->Add(this->label9);
			this->MainTask->Controls->Add(this->label10);
			this->MainTask->Controls->Add(this->textBox10);
			this->MainTask->Controls->Add(this->textBox11);
			this->MainTask->Controls->Add(this->label11);
			this->MainTask->Controls->Add(this->label12);
			this->MainTask->Controls->Add(this->textBox12);
			this->MainTask->Location = System::Drawing::Point(4, 22);
			this->MainTask->Name = L"MainTask";
			this->MainTask->Padding = System::Windows::Forms::Padding(3);
			this->MainTask->Size = System::Drawing::Size(1586, 653);
			this->MainTask->TabIndex = 1;
			this->MainTask->Text = L"Main Task";
			this->MainTask->UseVisualStyleBackColor = true;
			// 
			// textBox13
			// 
			this->textBox13->Location = System::Drawing::Point(126, 446);
			this->textBox13->Name = L"textBox13";
			this->textBox13->Size = System::Drawing::Size(48, 20);
			this->textBox13->TabIndex = 34;
			this->textBox13->Text = L"0";
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Location = System::Drawing::Point(88, 451);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(23, 13);
			this->label13->TabIndex = 33;
			this->label13->Text = L"U\'0";
			// 
			// dataGridView2
			// 
			this->dataGridView2->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView2->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(9) {
				this->dataGridViewTextBoxColumn1,
					this->dataGridViewTextBoxColumn2, this->dataGridViewTextBoxColumn3, this->dataGridViewTextBoxColumn4, this->dataGridViewTextBoxColumn5,
					this->dataGridViewTextBoxColumn6, this->dataGridViewTextBoxColumn7, this->dataGridViewTextBoxColumn8, this->dataGridViewTextBoxColumn9
			});
			this->dataGridView2->Location = System::Drawing::Point(637, 0);
			this->dataGridView2->Name = L"dataGridView2";
			this->dataGridView2->RowHeadersVisible = false;
			this->dataGridView2->Size = System::Drawing::Size(934, 396);
			this->dataGridView2->TabIndex = 32;
			this->dataGridView2->CellContentClick += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MyForm::dataGridView2_CellContentClick);
			// 
			// dataGridViewTextBoxColumn1
			// 
			this->dataGridViewTextBoxColumn1->HeaderText = L"Iter";
			this->dataGridViewTextBoxColumn1->Name = L"dataGridViewTextBoxColumn1";
			this->dataGridViewTextBoxColumn1->Width = 40;
			// 
			// dataGridViewTextBoxColumn2
			// 
			this->dataGridViewTextBoxColumn2->HeaderText = L"X";
			this->dataGridViewTextBoxColumn2->Name = L"dataGridViewTextBoxColumn2";
			this->dataGridViewTextBoxColumn2->ReadOnly = true;
			this->dataGridViewTextBoxColumn2->Width = 70;
			// 
			// dataGridViewTextBoxColumn3
			// 
			this->dataGridViewTextBoxColumn3->HeaderText = L"V_full";
			this->dataGridViewTextBoxColumn3->Name = L"dataGridViewTextBoxColumn3";
			this->dataGridViewTextBoxColumn3->ReadOnly = true;
			// 
			// dataGridViewTextBoxColumn4
			// 
			this->dataGridViewTextBoxColumn4->HeaderText = L"V_half2";
			this->dataGridViewTextBoxColumn4->Name = L"dataGridViewTextBoxColumn4";
			this->dataGridViewTextBoxColumn4->ReadOnly = true;
			// 
			// dataGridViewTextBoxColumn5
			// 
			this->dataGridViewTextBoxColumn5->HeaderText = L"V_full - V_half2";
			this->dataGridViewTextBoxColumn5->Name = L"dataGridViewTextBoxColumn5";
			// 
			// dataGridViewTextBoxColumn6
			// 
			this->dataGridViewTextBoxColumn6->HeaderText = L"ОЛП";
			this->dataGridViewTextBoxColumn6->Name = L"dataGridViewTextBoxColumn6";
			// 
			// dataGridViewTextBoxColumn7
			// 
			this->dataGridViewTextBoxColumn7->HeaderText = L"h";
			this->dataGridViewTextBoxColumn7->Name = L"dataGridViewTextBoxColumn7";
			this->dataGridViewTextBoxColumn7->Width = 70;
			// 
			// dataGridViewTextBoxColumn8
			// 
			this->dataGridViewTextBoxColumn8->HeaderText = L"C1";
			this->dataGridViewTextBoxColumn8->Name = L"dataGridViewTextBoxColumn8";
			this->dataGridViewTextBoxColumn8->Width = 40;
			// 
			// dataGridViewTextBoxColumn9
			// 
			this->dataGridViewTextBoxColumn9->HeaderText = L"C2";
			this->dataGridViewTextBoxColumn9->Name = L"dataGridViewTextBoxColumn9";
			this->dataGridViewTextBoxColumn9->Width = 40;
			// 
			// checkBox2
			// 
			this->checkBox2->AutoSize = true;
			this->checkBox2->Location = System::Drawing::Point(187, 446);
			this->checkBox2->Name = L"checkBox2";
			this->checkBox2->Size = System::Drawing::Size(76, 17);
			this->checkBox2->TabIndex = 31;
			this->checkBox2->Text = L"Fixed Step";
			this->checkBox2->UseVisualStyleBackColor = true;
			// 
			// button3
			// 
			this->button3->Location = System::Drawing::Point(532, 417);
			this->button3->Name = L"button3";
			this->button3->Size = System::Drawing::Size(80, 29);
			this->button3->TabIndex = 24;
			this->button3->Text = L"Zoom";
			this->button3->UseVisualStyleBackColor = true;
			this->button3->Click += gcnew System::EventHandler(this, &MyForm::button3_Click);
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(362, 425);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(34, 13);
			this->label7->TabIndex = 30;
			this->label7->Text = L"Nmax";
			// 
			// textBox7
			// 
			this->textBox7->Location = System::Drawing::Point(24, 420);
			this->textBox7->Name = L"textBox7";
			this->textBox7->Size = System::Drawing::Size(48, 20);
			this->textBox7->TabIndex = 19;
			this->textBox7->Text = L"0";
			// 
			// button4
			// 
			this->button4->Location = System::Drawing::Point(456, 417);
			this->button4->Name = L"button4";
			this->button4->Size = System::Drawing::Size(70, 29);
			this->button4->TabIndex = 17;
			this->button4->Text = L"Draw";
			this->button4->UseVisualStyleBackColor = true;
			this->button4->Click += gcnew System::EventHandler(this, &MyForm::button4_Click);
			// 
			// textBox8
			// 
			this->textBox8->Location = System::Drawing::Point(393, 421);
			this->textBox8->Name = L"textBox8";
			this->textBox8->Size = System::Drawing::Size(49, 20);
			this->textBox8->TabIndex = 29;
			this->textBox8->Text = L"250";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Location = System::Drawing::Point(5, 423);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(13, 13);
			this->label8->TabIndex = 18;
			this->label8->Text = L"a";
			// 
			// textBox9
			// 
			this->textBox9->Location = System::Drawing::Point(307, 421);
			this->textBox9->Name = L"textBox9";
			this->textBox9->Size = System::Drawing::Size(49, 20);
			this->textBox9->TabIndex = 28;
			this->textBox9->Text = L"0,01";
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Location = System::Drawing::Point(5, 447);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(13, 13);
			this->label9->TabIndex = 20;
			this->label9->Text = L"b";
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Location = System::Drawing::Point(277, 424);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(25, 13);
			this->label10->TabIndex = 27;
			this->label10->Text = L"Eps";
			// 
			// textBox10
			// 
			this->textBox10->Location = System::Drawing::Point(24, 444);
			this->textBox10->Name = L"textBox10";
			this->textBox10->Size = System::Drawing::Size(49, 20);
			this->textBox10->TabIndex = 21;
			this->textBox10->Text = L"2,0";
			// 
			// textBox11
			// 
			this->textBox11->Location = System::Drawing::Point(126, 420);
			this->textBox11->Name = L"textBox11";
			this->textBox11->Size = System::Drawing::Size(48, 20);
			this->textBox11->TabIndex = 26;
			this->textBox11->Text = L"10,0";
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Location = System::Drawing::Point(187, 423);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(13, 13);
			this->label11->TabIndex = 22;
			this->label11->Text = L"h";
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Location = System::Drawing::Point(88, 425);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(21, 13);
			this->label12->TabIndex = 25;
			this->label12->Text = L"U0";
			// 
			// textBox12
			// 
			this->textBox12->Location = System::Drawing::Point(206, 420);
			this->textBox12->Name = L"textBox12";
			this->textBox12->Size = System::Drawing::Size(61, 20);
			this->textBox12->TabIndex = 23;
			this->textBox12->Text = L"0,01";
			// 
			// zedGraphControl2
			// 
			this->zedGraphControl2->Location = System::Drawing::Point(0, 0);
			this->zedGraphControl2->Name = L"zedGraphControl2";
			this->zedGraphControl2->ScrollGrace = 0;
			this->zedGraphControl2->ScrollMaxX = 0;
			this->zedGraphControl2->ScrollMaxY = 0;
			this->zedGraphControl2->ScrollMaxY2 = 0;
			this->zedGraphControl2->ScrollMinX = 0;
			this->zedGraphControl2->ScrollMinY = 0;
			this->zedGraphControl2->ScrollMinY2 = 0;
			this->zedGraphControl2->Size = System::Drawing::Size(623, 370);
			this->zedGraphControl2->TabIndex = 1;
			this->zedGraphControl2->Load += gcnew System::EventHandler(this, &MyForm::zedGraphControl2_Load);
			// 
			// tabControl2
			// 
			this->tabControl2->Controls->Add(this->ux);
			this->tabControl2->Controls->Add(this->uproizvx);
			this->tabControl2->Controls->Add(this->phase);
			this->tabControl2->Location = System::Drawing::Point(0, 0);
			this->tabControl2->Name = L"tabControl2";
			this->tabControl2->SelectedIndex = 0;
			this->tabControl2->Size = System::Drawing::Size(631, 396);
			this->tabControl2->TabIndex = 35;
			// 
			// ux
			// 
			this->ux->Controls->Add(this->zedGraphControl2);
			this->ux->Location = System::Drawing::Point(4, 22);
			this->ux->Name = L"ux";
			this->ux->Padding = System::Windows::Forms::Padding(3);
			this->ux->Size = System::Drawing::Size(623, 370);
			this->ux->TabIndex = 0;
			this->ux->Text = L"U(x)";
			this->ux->UseVisualStyleBackColor = true;
			// 
			// uproizvx
			// 
			this->uproizvx->Controls->Add(this->zedGraphControl3);
			this->uproizvx->Location = System::Drawing::Point(4, 22);
			this->uproizvx->Name = L"uproizvx";
			this->uproizvx->Padding = System::Windows::Forms::Padding(3);
			this->uproizvx->Size = System::Drawing::Size(623, 370);
			this->uproizvx->TabIndex = 1;
			this->uproizvx->Text = L"U\'(x)";
			this->uproizvx->UseVisualStyleBackColor = true;
			// 
			// phase
			// 
			this->phase->Controls->Add(this->zedGraphControl4);
			this->phase->Location = System::Drawing::Point(4, 22);
			this->phase->Name = L"phase";
			this->phase->Size = System::Drawing::Size(623, 370);
			this->phase->TabIndex = 2;
			this->phase->Text = L"Фазовый Портрет";
			this->phase->UseVisualStyleBackColor = true;
			// 
			// zedGraphControl3
			// 
			this->zedGraphControl3->Location = System::Drawing::Point(0, 0);
			this->zedGraphControl3->Name = L"zedGraphControl3";
			this->zedGraphControl3->ScrollGrace = 0;
			this->zedGraphControl3->ScrollMaxX = 0;
			this->zedGraphControl3->ScrollMaxY = 0;
			this->zedGraphControl3->ScrollMaxY2 = 0;
			this->zedGraphControl3->ScrollMinX = 0;
			this->zedGraphControl3->ScrollMinY = 0;
			this->zedGraphControl3->ScrollMinY2 = 0;
			this->zedGraphControl3->Size = System::Drawing::Size(623, 370);
			this->zedGraphControl3->TabIndex = 36;
			// 
			// zedGraphControl4
			// 
			this->zedGraphControl4->Location = System::Drawing::Point(0, 0);
			this->zedGraphControl4->Name = L"zedGraphControl4";
			this->zedGraphControl4->ScrollGrace = 0;
			this->zedGraphControl4->ScrollMaxX = 0;
			this->zedGraphControl4->ScrollMaxY = 0;
			this->zedGraphControl4->ScrollMaxY2 = 0;
			this->zedGraphControl4->ScrollMinX = 0;
			this->zedGraphControl4->ScrollMinY = 0;
			this->zedGraphControl4->ScrollMinY2 = 0;
			this->zedGraphControl4->Size = System::Drawing::Size(623, 370);
			this->zedGraphControl4->TabIndex = 37;
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1606, 680);
			this->Controls->Add(this->tabControl1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->tabControl1->ResumeLayout(false);
			this->TestTask->ResumeLayout(false);
			this->TestTask->PerformLayout();
			this->MainTask->ResumeLayout(false);
			this->MainTask->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView2))->EndInit();
			this->tabControl2->ResumeLayout(false);
			this->ux->ResumeLayout(false);
			this->uproizvx->ResumeLayout(false);
			this->phase->ResumeLayout(false);
			this->ResumeLayout(false);

		}
#pragma endregion
	private: 
		double f1(double x){
			return sin(x);
		}

		double f2(double x) {
			return sin(2 * x);
		}

	// Test task 
	private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
		/* --- Нажали Draw (button 1) --- */

		GraphPane^ panel = zedGraphControl1->GraphPane;
		panel->CurveList->Clear();
		PointPairList^ numeric_trajectory_list = gcnew ZedGraph::PointPairList();
		PointPairList^ real_trajectory_list = gcnew ZedGraph::PointPairList();

		// Считываем данные, просчитываем то, что нужно double x0, double U0, double h, double b, int Nmax

		double xmin = Convert::ToDouble(textBox1->Text);
		double xmax = Convert::ToDouble(textBox2->Text);
		double U0 = Convert::ToDouble(textBox5->Text);
		double h = Convert::ToDouble(textBox3->Text);
		double Eps = Convert::ToDouble(textBox4->Text);
		double Nmax = Convert::ToDouble(textBox6->Text);

		std::vector<DataTestTask> TestData;

		bool isFixedStep = checkBox1->Checked;

		// Просчитываем численную траекторию методом RK4
		if (isFixedStep) {
			TestData = RK4_method_fixed_step_test_task(xmin, U0, h, xmax, Nmax);
		}
		else {
			TestData = RK4_method_addaptive_step_test_task(xmin, U0, h, xmax, Nmax, Eps);
		}


		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;
		
		// Вывол данных
		dataGridView1->Rows->Clear();
		
		for (int i = 0; i < TestData.size(); i++)
		{
			DataTestTask currentStep = TestData[i];
			//Добавление на график точек численной траектории, реальной траектории
			numeric_trajectory_list->Add(currentStep.x, currentStep.V_full);
			real_trajectory_list->Add(currentStep.x, currentStep.U);

			//Печать в таблицу
			DataGridViewCellStyle^ style = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			style->Format = L"F10";

			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i;

			dataGridView1->Rows[i]->Cells[1]->Value = currentStep.x; 	
			//dataGridView1->Rows[i]->Cells[0]->Style = style;

			dataGridView1->Rows[i]->Cells[2]->Value = currentStep.V_full;
			dataGridView1->Rows[i]->Cells[2]->Style = style;

			dataGridView1->Rows[i]->Cells[3]->Value = currentStep.V_half2;
			dataGridView1->Rows[i]->Cells[3]->Style = style;

			dataGridView1->Rows[i]->Cells[4]->Value = currentStep.Vfull_Vhalf2;
			dataGridView1->Rows[i]->Cells[4]->Style = style;

			dataGridView1->Rows[i]->Cells[5]->Value = currentStep.OLP;
			dataGridView1->Rows[i]->Cells[5]->Style = style;

			dataGridView1->Rows[i]->Cells[6]->Value = currentStep.h;
			dataGridView1->Rows[i]->Cells[6]->Style = style;

			dataGridView1->Rows[i]->Cells[7]->Value = currentStep.C1;
			//dataGridView1->Rows[i]->Cells[0]->Style = style;

			dataGridView1->Rows[i]->Cells[8]->Value = currentStep.C2;
			// dataGridView1->Rows[i]->Cells[0]->Style = style;

			dataGridView1->Rows[i]->Cells[9]->Value = currentStep.U;
			dataGridView1->Rows[i]->Cells[9]->Style = style;

			dataGridView1->Rows[i]->Cells[10]->Value = currentStep.U_V;
			dataGridView1->Rows[i]->Cells[10]->Style = style;


		}
		LineItem Curve1 = panel->AddCurve("V", numeric_trajectory_list, Color::Red,SymbolType::Plus);
		LineItem Curve2 = panel->AddCurve("U", real_trajectory_list, Color::Blue, SymbolType::None);

		// Устанавливаем интересующий нас интервал по оси X
		panel->XAxis->Scale->Min = xmin_limit;
		panel->XAxis->Scale->Max = xmax_limit;
/*
		// Устанавливаем интересующий нас интервал по оси Y
		panel->YAxis->Scale->Min = ymin_limit;
		panel->YAxis->Scale->Max = ymax_limit;
*/
		// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
		// В противном случае на рисунке будет показана только часть графика, 
		// которая умещается в интервалы по осям, установленные по умолчанию
		zedGraphControl1->AxisChange();
		// Обновляем график
		zedGraphControl1->Invalidate();

	}

	// Main task
	private: System::Void button4_Click(System::Object^ sender, System::EventArgs^ e) {
		/* --- Нажали Draw (button 4) --- */

		GraphPane^ panel_ux = zedGraphControl2->GraphPane;
		GraphPane^ panel_uproizvx = zedGraphControl3->GraphPane;
		GraphPane^ panel_phase = zedGraphControl4->GraphPane;


		panel_ux->CurveList->Clear();
		panel_uproizvx->CurveList->Clear();
		panel_phase->CurveList->Clear();


		PointPairList^ numeric_trajectory_list = gcnew ZedGraph::PointPairList();
		PointPairList^ proizv_trajectory_list = gcnew ZedGraph::PointPairList();
		PointPairList^ phase_list = gcnew ZedGraph::PointPairList();



		// Считываем данные, просчитываем то, что нужно double x0, double U1_0, double U2_0, double h, double b, int Nmax, double Eps

		double xmin = Convert::ToDouble(textBox7->Text);
		double xmax = Convert::ToDouble(textBox10->Text);
		double U1_0 = Convert::ToDouble(textBox11->Text);
		double U2_0 = Convert::ToDouble(textBox13->Text);

		double h = Convert::ToDouble(textBox12->Text);
		double Eps = Convert::ToDouble(textBox9->Text);
		double Nmax = Convert::ToDouble(textBox8->Text);

		std::vector<DataMainTask> TestData;

		bool isFixedStep = checkBox2->Checked;

		// Просчитываем численную траекторию методом RK4
		if (isFixedStep) {
			TestData = RK4_method_fixed_step_main_task(xmin, U1_0, U2_0, h, xmax, Nmax);
		}
		else {
			TestData = RK4_method_addaptive_step_main_task(xmin, U1_0, U2_0, h, xmax, Nmax, Eps);
		}


		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;

		// Вывол данных
		dataGridView2->Rows->Clear();

		for (int i = 0; i < TestData.size(); i++)
		{
			DataMainTask currentStep = TestData[i];
			//Добавление на график точек численной траектории
			numeric_trajectory_list->Add(currentStep.x, currentStep.V_full);
			proizv_trajectory_list->Add(currentStep.x, currentStep.proizv);
			phase_list->Add(currentStep.V_full, currentStep.proizv);



			//Печать в таблицу
			DataGridViewCellStyle^ style = (gcnew System::Windows::Forms::DataGridViewCellStyle());
			style->Format = L"F10";

			dataGridView2->Rows->Add();
			dataGridView2->Rows[i]->Cells[0]->Value = i;

			dataGridView2->Rows[i]->Cells[1]->Value = currentStep.x;
			//dataGridView2->Rows[i]->Cells[0]->Style = style;

			dataGridView2->Rows[i]->Cells[2]->Value = currentStep.V_full;
			dataGridView2->Rows[i]->Cells[2]->Style = style;

			dataGridView2->Rows[i]->Cells[3]->Value = currentStep.V_half2;
			dataGridView2->Rows[i]->Cells[3]->Style = style;

			dataGridView2->Rows[i]->Cells[4]->Value = currentStep.Vfull_Vhalf2;
			dataGridView2->Rows[i]->Cells[4]->Style = style;

			dataGridView2->Rows[i]->Cells[5]->Value = currentStep.OLP;
			dataGridView2->Rows[i]->Cells[5]->Style = style;

			dataGridView2->Rows[i]->Cells[6]->Value = currentStep.h;
			dataGridView2->Rows[i]->Cells[6]->Style = style;

			dataGridView2->Rows[i]->Cells[7]->Value = currentStep.C1;
			//dataGridView2->Rows[i]->Cells[0]->Style = style;

			dataGridView2->Rows[i]->Cells[8]->Value = currentStep.C2;
			// dataGridView2->Rows[i]->Cells[0]->Style = style;



		}
		LineItem Curve1 = panel_ux->AddCurve("U(x)", numeric_trajectory_list, Color::Red, SymbolType::Plus);
		LineItem Curve2 = panel_uproizvx->AddCurve("U'(x)", proizv_trajectory_list, Color::Red, SymbolType::Plus);
		LineItem Curvee = panel_phase->AddCurve("U' vs U", phase_list, Color::Red, SymbolType::Plus);



		// Устанавливаем интересующий нас интервал по оси X
		panel_ux->XAxis->Scale->Min = xmin_limit;
		panel_ux->XAxis->Scale->Max = xmax_limit;
		
		// Автоматический масштаб 
		panel_uproizvx = zedGraphControl3->GraphPane;
		panel_uproizvx->YAxis->Scale->MinAuto = true;
		panel_uproizvx->YAxis->Scale->MaxAuto = true;
		panel_uproizvx->XAxis->Scale->MinAuto = true;
		panel_uproizvx->XAxis->Scale->MaxAuto = true;

		// Подписи осей
		panel_uproizvx->XAxis->Title->Text = "x";
		panel_uproizvx->YAxis->Title->Text = "U'(x)";
		panel_uproizvx->Title->Text = "Velocity vs Time";

		panel_phase = zedGraphControl4->GraphPane;

		// Авто-масштаб для обеих осей
		panel_phase->XAxis->Scale->MinAuto = true;
		panel_phase->XAxis->Scale->MaxAuto = true;
		panel_phase->YAxis->Scale->MinAuto = true;
		panel_phase->YAxis->Scale->MaxAuto = true;

		// Подписи осей
		panel_phase->XAxis->Title->Text = "U (displacement)";
		panel_phase->YAxis->Title->Text = "U' (velocity)";
		panel_phase->Title->Text = "Phase Portrait";

		panel_phase->XAxis->Scale->MajorStepAuto = true;  // автоматический шаг сетки
		panel_phase->YAxis->Scale->MajorStepAuto = true;
		panel_phase->XAxis->Scale->MinorStepAuto = true;
		panel_phase->YAxis->Scale->MinorStepAuto = true;
		/*
				// Устанавливаем интересующий нас интервал по оси Y
				panel->YAxis->Scale->Min = ymin_limit;
				panel->YAxis->Scale->Max = ymax_limit;
		*/
		// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
		// В противном случае на рисунке будет показана только часть графика, 
		// которая умещается в интервалы по осям, установленные по умолчанию
		zedGraphControl2->AxisChange();
		// Обновляем график
		zedGraphControl2->Invalidate();

	}

	private: System::Void zedGraphControl1_Load(System::Object^  sender, System::EventArgs^  e) {
	}

private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
	
	GraphPane^ panel = zedGraphControl1->GraphPane;
	double xmin = Convert::ToDouble(textBox5->Text);
	double xmax = Convert::ToDouble(textBox4->Text);
	// Устанавливаем интересующий нас интервал по оси X
	panel->XAxis->Scale->Min = xmin;
	panel->XAxis->Scale->Max = xmax;

	// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
	// В противном случае на рисунке будет показана только часть графика, 
	// которая умещается в интервалы по осям, установленные по умолчанию
	zedGraphControl1->AxisChange();
	// Обновляем график
	zedGraphControl1->Invalidate();

}
private: System::Void button3_Click(System::Object^ sender, System::EventArgs^ e) {

	GraphPane^ panel = zedGraphControl2->GraphPane;
	double xmin = Convert::ToDouble(textBox11->Text);
	double xmax = Convert::ToDouble(textBox9->Text);
	// Устанавливаем интересующий нас интервал по оси X
	panel->XAxis->Scale->Min = xmin;
	panel->XAxis->Scale->Max = xmax;

	// Вызываем метод AxisChange (), чтобы обновить данные об осях. 
	// В противном случае на рисунке будет показана только часть графика, 
	// которая умещается в интервалы по осям, установленные по умолчанию
	zedGraphControl2->AxisChange();
	// Обновляем график
	zedGraphControl2->Invalidate();

}
private: System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void dataGridView1_CellContentClick(System::Object^ sender, System::Windows::Forms::DataGridViewCellEventArgs^ e) {
}
private: System::Void textBox1_TextChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label1_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label5_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label4_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void label6_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void checkBox1_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void TestTask_Click(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void zedGraphControl2_Load(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void dataGridView2_CellContentClick(System::Object^ sender, System::Windows::Forms::DataGridViewCellEventArgs^ e) {
}
};
}
