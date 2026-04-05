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
double variant = 1;
double f_test_task(double x, double u) { 
	double l = ((int)variant % 2 ? -1 : 1) * (variant / 2.0);
	return l * u;
}
static double f_test_pervoobr(double x, double u) {
	double lambda = ((static_cast<int>(variant) % 2 == 0) ? 1.0 : -1.0) * (variant / 2.0);
	return u * exp(lambda * x);
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

		V_n = RK4_step_test_task(V_n, x_n, h);
		x_n += h;

		// Подсчет двойного шага для таблицы
		double V_half = RK4_step_test_task(V_n, x_n, h / 2.0);
		double V_half2 = RK4_step_test_task(V_half, x_n + h / 2.0, h / 2.0);

		CurrStepData.x = x_n;
		CurrStepData.V_full = V_n;
		CurrStepData.V_half2 = V_half2;
		CurrStepData.Vfull_Vhalf2 = V_n - V_half2;
		CurrStepData.OLP = (V_n - V_half2) / 15.0;
		CurrStepData.h = h;
		CurrStepData.C1 = 0;
		CurrStepData.C2 = 0;
		CurrStepData.U = f_test_pervoobr(x_n, V_n);
		CurrStepData.U = f_test_pervoobr(x_n, V_n) - V_n;

		Data.push_back(CurrStepData);
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
		CurrStepData.OLP = (V_full - V_half2) / 15.0;
		CurrStepData.h = h_n;
		CurrStepData.C1 = 0;
		CurrStepData.C2 = 0;
		CurrStepData.U = f_test_pervoobr(x_n, V_n);
		CurrStepData.U = f_test_pervoobr(x_n, V_n) - V_n;

		Data.push_back(CurrStepData);
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

		V_n = V_n1;
		x_n += h;

		Data.push_back(CurrStepData);
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
		CurrStepData.C1 = 0;
		CurrStepData.C2 = 0;

		Data.push_back(CurrStepData);
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
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  X;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  F_1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^  F_2;
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
			this->X = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->F_1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->F_2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
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
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->SuspendLayout();
			// 
			// zedGraphControl1
			// 
			this->zedGraphControl1->Location = System::Drawing::Point(38, 30);
			this->zedGraphControl1->Name = L"zedGraphControl1";
			this->zedGraphControl1->ScrollGrace = 0;
			this->zedGraphControl1->ScrollMaxX = 0;
			this->zedGraphControl1->ScrollMaxY = 0;
			this->zedGraphControl1->ScrollMaxY2 = 0;
			this->zedGraphControl1->ScrollMinX = 0;
			this->zedGraphControl1->ScrollMinY = 0;
			this->zedGraphControl1->ScrollMinY2 = 0;
			this->zedGraphControl1->Size = System::Drawing::Size(501, 327);
			this->zedGraphControl1->TabIndex = 0;
			this->zedGraphControl1->Load += gcnew System::EventHandler(this, &MyForm::zedGraphControl1_Load);
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(633, 386);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(142, 29);
			this->button1->TabIndex = 1;
			this->button1->Text = L"Draw";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// dataGridView1
			// 
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(3) {
				this->X, this->F_1,
					this->F_2
			});
			this->dataGridView1->Location = System::Drawing::Point(559, 30);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->RowHeadersVisible = false;
			this->dataGridView1->Size = System::Drawing::Size(274, 327);
			this->dataGridView1->TabIndex = 2;
			// 
			// X
			// 
			this->X->HeaderText = L"X";
			this->X->Name = L"X";
			this->X->ReadOnly = true;
			this->X->Width = 50;
			// 
			// F_1
			// 
			this->F_1->HeaderText = L"F_1";
			this->F_1->Name = L"F_1";
			this->F_1->ReadOnly = true;
			// 
			// F_2
			// 
			this->F_2->HeaderText = L"F_2";
			this->F_2->Name = L"F_2";
			this->F_2->ReadOnly = true;
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(59, 394);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(13, 13);
			this->label1->TabIndex = 3;
			this->label1->Text = L"a";
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(78, 394);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(48, 20);
			this->textBox1->TabIndex = 4;
			this->textBox1->Text = L"0";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(171, 396);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(13, 13);
			this->label2->TabIndex = 5;
			this->label2->Text = L"b";
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(190, 393);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(49, 20);
			this->textBox2->TabIndex = 6;
			this->textBox2->Text = L"1";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(287, 398);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(13, 13);
			this->label3->TabIndex = 7;
			this->label3->Text = L"h";
			// 
			// textBox3
			// 
			this->textBox3->Location = System::Drawing::Point(306, 394);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(61, 20);
			this->textBox3->TabIndex = 8;
			this->textBox3->Text = L"0,1";
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(633, 437);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(142, 29);
			this->button2->TabIndex = 9;
			this->button2->Text = L"Zoom";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// textBox4
			// 
			this->textBox4->Location = System::Drawing::Point(190, 437);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(49, 20);
			this->textBox4->TabIndex = 13;
			this->textBox4->Text = L"1";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(171, 440);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(13, 13);
			this->label4->TabIndex = 12;
			this->label4->Text = L"b";
			// 
			// textBox5
			// 
			this->textBox5->Location = System::Drawing::Point(78, 436);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(48, 20);
			this->textBox5->TabIndex = 11;
			this->textBox5->Text = L"0";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(59, 438);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(13, 13);
			this->label5->TabIndex = 10;
			this->label5->Text = L"a";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(858, 497);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->textBox5);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->textBox3);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->Controls->Add(this->dataGridView1);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->zedGraphControl1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: 
		double f1(double x){
			return sin(x);
		}

		double f2(double x) {
			return sin(2 * x);
		}

	private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {

		GraphPane^ panel = zedGraphControl1->GraphPane;
		panel->CurveList->Clear();
		PointPairList^ f1_list = gcnew ZedGraph::PointPairList();
		PointPairList^ f2_list = gcnew ZedGraph::PointPairList();

		// Интервал, где есть данные
		double xmin = Convert::ToDouble(textBox1->Text);
		double xmax = Convert::ToDouble(textBox2->Text);

		double h = Convert::ToDouble(textBox3->Text);


		double xmin_limit = xmin - 0.1;
		double xmax_limit = xmax + 0.1;
/*
		double ymin_limit = -1.0;
		double ymax_limit = 100.0;
*/
		// Список точек
		int i = 0;
		dataGridView1->Rows->Clear();
		for (double x = xmin; x <= xmax; x += h)
		{
			//Добавление на график
			f1_list->Add(x, f1(x));
			f2_list->Add(x, f2(x));
			//Печать в таблицу
			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = x; 			
			dataGridView1->Rows[i]->Cells[1]->Value = floor(f1(x) * 1000) / 1000;
			dataGridView1->Rows[i]->Cells[2]->Value = floor(f2(x) * 1000) / 1000;
			i++;
		}
		LineItem Curve1 = panel->AddCurve("F1(x)", f1_list, Color::Red,SymbolType::Plus);
		LineItem Curve2 = panel->AddCurve("F2(x)", f2_list, Color::Blue, SymbolType::None);

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
private: System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
}
};
}
