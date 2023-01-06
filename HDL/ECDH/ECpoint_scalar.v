module ECpoint_scalar (
	input clk,
	input rst_n,
	input start,
	input [`BW_GF-1:0] k,
	input [`BW_GF-1:0] Px,
	input [`BW_GF-1:0] Py,
	output [`BW_GF-1:0] Qx,
	output [`BW_GF-1:0] Qy,
	output reg valid
);

integer i;

// state
localparam IDLE = 0, RESET = 1, SEARCH = 2, INIT = 3;
localparam PT_DBL = 4, PT_ADD = 5, CALXSQ = 6, MULINV = 7;

// Register Alias
localparam A = 0, B = 1, C = 2, D = 3, E = 4, F = 5, G = 6, H = 7, I = 8;

// State
reg valid_n;
reg [2:0] state, state_n;
reg is_first, is_first_n;
reg [$clog2(`BW_GF)-1:0] cnt, cnt_n;
reg done_pointcal, done_pointcal_n;
reg kbit, kbit_n;
reg [3:0] step_cnt, step_cnt_n;
reg ctrl_cnt, ctrl_cnt_n;
reg done_calc;	// also serve as write enable

// data register: A~I
reg [`BW_GF-1:0] data[0:9-1];
reg [`BW_GF-1:0] data_n[0:9-1];

// 	control signal
reg [4-1:0] sel_mul_A, sel_mul_A_n;
reg [3-1:0] sel_mul_B, sel_mul_B_n;
reg [3-1:0] sel_add_A, sel_add_A_n;
reg [3-1:0] sel_add_B, sel_add_B_n;
reg [2-1:0] w_selA, w_selA_tmp;
reg [2-1:0] w_selB, w_selB_tmp;
reg [2-1:0] w_selC, w_selC_tmp;
reg [2-1:0] w_selD, w_selD_tmp;
reg [2-1:0] w_selE, w_selE_tmp;
reg [2-1:0] w_selF, w_selF_tmp;
reg [2-1:0] w_selG, w_selG_tmp;
reg [1-1:0] w_selH, w_selH_tmp;
reg [2-1:0] w_selI, w_selI_tmp;
reg en_mul, en_mul_n;
reg en_add, en_add_n;
reg en_multinv, en_multinv_n;
reg is_sub, is_sub_n;

// data flow
reg [`BW_GF-1:0] mult1, mult2;
reg [`BW_GF-1:0] add1, add2;
wire [`BW_GF-1:0] sum;
wire [`BW_GF-1:0] prod;
wire [`BW_GF-1:0] z_multinv;
wire [9-1:0] powerk;
reg [3:0] power_cnt, power_cnt_n;
reg inv2power_mode, inv2power_mode_n;	// 0: double, 1: addition
// reg done_inv2power;
wire [10-1:0] expshift;

// valid signal
wire valid_mod_prod;
wire valid_mod_sum;
wire valid_multinv;

assign Qx = data[A];
assign Qy = data[B];

always @(*) 
	if (state == MULINV & done_pointcal)
		valid_n = 1;
	else
		valid_n = 0;

always @(posedge clk) begin
	valid <= valid_n;
end

// data register file ---------------------------------------------------------
always @(posedge clk)
	for (i = 0; i < 9; i = i + 1)
		data[i] <= data_n[i];

always @(*) begin
	for (i = 0; i < 9; i = i + 1)
		data_n[i] = data[i];
	case(w_selA)
		2'h1:	data_n[A] = prod;
		2'h2:	data_n[A] = sum;
		2'h3:	data_n[A] = Px;
	endcase
	case(w_selB)
		2'h1:	data_n[B] = prod;
		2'h2:	data_n[B] = sum;
		2'h3:	data_n[B] = Py;
	endcase
	case(w_selC)
		2'h1:	data_n[C] = prod;
		2'h2:	data_n[C] = sum;
		2'h3:	data_n[C] = 1;
	endcase
	case(w_selD)
		2'h1:	data_n[D] = prod;
		2'h2:	data_n[D] = sum;
		2'h3:	data_n[D] = z_multinv;
	endcase
	case(w_selE)
		2'h1:	data_n[E] = prod;
		2'h2:	data_n[E] = sum;
	endcase
	case(w_selF)
		2'h1:	data_n[F] = prod;
		2'h2:	data_n[F] = sum;
	endcase
	case(w_selG)
		2'h1:	data_n[G] = prod;
		2'h3:	data_n[G] = Px;
	endcase
	case(w_selH)
		1'h1:	data_n[H] = Py;
	endcase
	case(w_selI)
		2'h2:	data_n[I] = sum;
		2'h3:	data_n[I] = 1;
	endcase
end

// Multiplication and Addition Input selector ---------------------------------
always @(*) begin
	case(sel_mul_A)
		4'h0:	mult1 = data[A];
		4'h1:	mult1 = data[B];
		4'h2:	mult1 = data[C];
		4'h3:	mult1 = data[D];
		4'h4:	mult1 = data[E];
		4'h5:	mult1 = data[F];
		4'h6:	mult1 = data[G];
		4'h7:	mult1 = data[H];
		4'h8:	mult1 = data[I];
		4'h9:	mult1 = 3;
		4'hA:	mult1 = 4;
		4'hB:	mult1 = 8;
		4'hC:	mult1 = `COEF_A;
		4'hD:	mult1 = `TWO_MI;
		default: mult1 = 0;
	endcase
	case(sel_mul_B)
		3'h0:	mult2 = data[A];
		3'h1:	mult2 = data[B];
		3'h2:	mult2 = data[C];
		3'h3:	mult2 = data[D];
		3'h4:	mult2 = data[E];
		3'h5:	mult2 = data[F];
		3'h6:	mult2 = data[I];
		default: mult2 = 0;
	endcase
	case(sel_add_A)
		3'h0:	add1 = data[A];
		3'h1:	add1 = data[B];
		3'h2:	add1 = data[C];
		3'h3:	add1 = data[D];
		3'h4:	add1 = data[E];
		3'h5:	add1 = data[F];
		3'h6:	add1 = data[G];
		3'h7:	add1 = `PRIME;
		default: add1 = 0;
	endcase
	// add1 = data[{1'b0, sel_add_A}];
	add2 = data[{1'b0, sel_add_B}];
end

// Multiplication and Module block ------------------------------------------
Multiplication_192x192 u_MULT(
	.clk(clk),
	.rst(en_mul),
	.a(mult1),
	.b(mult2),
	.out(prod),
	.valid(valid_mod_prod)
);

// Summation and Modulo block -----------------------------------------------
ADD_192 u_ADD(
	.clk(clk),
	.en(en_add),
	.a(add1),
	.b(add2),
	.is_sub(is_sub),
	.out(sum),
	.valid(valid_mod_sum)
);

// FSM ------------------------------------------------------------------------
always @(*) begin
	state_n = state;
	case(state)
		IDLE:	state_n = (start)? RESET: IDLE;
		RESET:	state_n = SEARCH;
		SEARCH:	// check scalar bit 
			if (cnt < `BW_GF) begin
				if (is_first == 1) begin
					if (kbit_n == 1) 
						state_n = INIT;
				end else
					state_n = PT_DBL;
			end else
				state_n = MULINV;
		INIT:	state_n = CALXSQ;
		CALXSQ:	state_n = (done_calc)? SEARCH: CALXSQ;
		PT_DBL:	
			if (done_pointcal)
				if (kbit)
					state_n = PT_ADD;
				else
					state_n = CALXSQ;
		PT_ADD:	state_n = (done_pointcal)? SEARCH: PT_ADD;
		MULINV:	state_n = (done_pointcal)? IDLE: MULINV;
	endcase
end

always @(posedge clk)
	if (~rst_n)
		state <= IDLE;
	else 
		state <= state_n;

always @(*) begin
	cnt_n = cnt;
	if (state == SEARCH)
		cnt_n = cnt + 1;
	else if (state == RESET)
		cnt_n = 0;
end

reg [`BW_GF-1:0] kshift;

always @(*) begin
	kbit_n = kbit;
	if (state == SEARCH) begin
		kshift = k << cnt;
		kbit_n = kshift[`BW_GF-1];
	end
end

always @(*) begin
	if (state == RESET)	is_first_n = 1;
	else if (state == INIT) is_first_n = 0;
	else is_first_n = is_first;
end

always @(posedge clk) begin
	cnt <= cnt_n;
	kbit <= kbit_n;
	is_first <= is_first_n;
end

// Point Computing Steps ------------------------------------------------------
// Inverse 2 exponential calculating
assign expshift = powerk << step_cnt;

always @(*) begin
	step_cnt_n = 0;
	done_pointcal_n = 0;
	ctrl_cnt_n = 0;
	inv2power_mode_n = inv2power_mode;
	case(state)
		PT_DBL: begin
			if (step_cnt < 13) begin
				if (done_calc) begin
					ctrl_cnt_n = 0;
					step_cnt_n = step_cnt + 1;
				end else begin
					ctrl_cnt_n = 1;
					step_cnt_n = step_cnt;	
				end
			end else begin
				step_cnt_n = 0;
				ctrl_cnt_n = 0;
			end
			if (step_cnt == 12 & done_calc)
				done_pointcal_n = 1;
		end
		PT_ADD:
			if (done_calc) begin
				if (step_cnt < 13) begin
					ctrl_cnt_n = 0;
					step_cnt_n = step_cnt + 1;	
				end else begin
					ctrl_cnt_n = 1;
					step_cnt_n = step_cnt;
					done_pointcal_n = 1;
				end
			end else begin
				ctrl_cnt_n = 1;
				step_cnt_n = step_cnt;
			end
		CALXSQ:
			ctrl_cnt_n = 1;
		MULINV: begin
			if (done_calc) begin	// if done one mission
				ctrl_cnt_n = 0;
				if (step_cnt == 0) begin				// AMI running
					step_cnt_n = step_cnt + 1;
					inv2power_mode_n = 1;
				end else if (step_cnt < 10) begin 		// calculate 2^(-k)
					if (inv2power_mode == 0 & expshift[9] == 1) begin
						inv2power_mode_n = 1;
						step_cnt_n = step_cnt;
					end else begin
						inv2power_mode_n = 0;
						step_cnt_n = step_cnt + 1;
					end
				end else if (step_cnt < 15) begin
					step_cnt_n = step_cnt + 1;
				end else begin
					ctrl_cnt_n = 1;
					done_pointcal_n = 1;
					step_cnt_n = step_cnt;
				end
			end else begin
				ctrl_cnt_n = 1;
				step_cnt_n = step_cnt;
			end
		end
		default: begin
			step_cnt_n = 0;
			done_pointcal_n = 0;	
			ctrl_cnt_n = 0;
			inv2power_mode_n = 0;
		end
	endcase
end

always @(posedge clk) begin
	step_cnt <= step_cnt_n;
	done_pointcal <= done_pointcal_n;
	ctrl_cnt <= ctrl_cnt_n;
	inv2power_mode <= inv2power_mode_n;
end

// Arithmetic Enable control --------------------------------------------------
always @(*) begin
	en_add_n = 0;
	en_mul_n = 0;
	case (state)
		PT_DBL: begin
			if (step_cnt == 0 & ctrl_cnt == 0)
				en_mul_n = 1;
			else if (step_cnt < 12 & valid_mod_prod)
				en_mul_n = 1;
			
			case(step_cnt)
				4, 5, 8, 9, 10, 12:
					en_add_n = ~ctrl_cnt;
			endcase
		end
		PT_ADD: begin
			if (done_calc) begin
				if (step_cnt == 3) 
					en_mul_n = 0;
				else 
					en_mul_n = 1;
			end

			case(step_cnt)
				3, 4, 6, 7, 8, 9, 10, 12:
					if (ctrl_cnt == 0)
						en_add_n = 1;
			endcase
		end
		CALXSQ: begin
			en_mul_n = ~ctrl_cnt;
		end
		MULINV: begin
			case(step_cnt)
				1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15:
					if (ctrl_cnt == 0) en_mul_n = 1;
				11:
					if (ctrl_cnt == 0) en_add_n = 1;
			endcase
		end
		default:	begin
			en_add_n = 0;
			en_mul_n = 0;
		end
	endcase
end

always @(*) begin
	if (state == MULINV & step_cnt == 0)
		en_multinv_n = ~ctrl_cnt;
	else
		en_multinv_n = 0;
end

always @(posedge clk) begin
	en_add <= en_add_n;
	en_mul <= en_mul_n;
	en_multinv <= en_multinv_n;
end

// Computing I/O control ------------------------------------------------------
always @(*) begin
	sel_add_A_n = 0;
	sel_add_B_n = 0;
	sel_mul_A_n = 0;
	sel_mul_B_n = 0;
	is_sub_n = 0;
	w_selA_tmp = 0;
    w_selB_tmp = 0;
    w_selC_tmp = 0;
    w_selD_tmp = 0;
    w_selE_tmp = 0;
    w_selF_tmp = 0;
    w_selG_tmp = 0;
    w_selH_tmp = 0;
    w_selI_tmp = 0;
	case(state)
		INIT:	begin
			w_selA_tmp = 3;
			w_selB_tmp = 3;
			w_selC_tmp = 3;
		end
		PT_DBL: 
			case(step_cnt)
				4'h0:	begin
					sel_mul_A_n = 9;
					sel_mul_B_n = 3;
					w_selE_tmp = 1;
				end
				4'h1:	begin
					sel_mul_A_n = 2;
					sel_mul_B_n = 2;
					w_selD_tmp = 1;
				end
				4'h2:	begin
					sel_mul_A_n = 3;
					sel_mul_B_n = 3;
					w_selD_tmp = 1;
				end
				4'h3:	begin
					sel_mul_A_n = 4'hC;
					sel_mul_B_n = 3;
					w_selD_tmp = 1;
				end
				4'h4:	begin
					sel_mul_A_n = 1;
					sel_mul_B_n = 2;
					w_selD_tmp = 1;
					sel_add_A_n = 3;
					sel_add_B_n = 4;
					w_selE_tmp = 2;
				end
				4'h5:	begin
					sel_mul_A_n = 1;
					sel_mul_B_n = 1;
					w_selB_tmp = 1;
					sel_add_A_n = 3;
					sel_add_B_n = 3;
					w_selC_tmp = 2;
				end
				4'h6:	begin
					sel_mul_A_n = 1;
					sel_mul_B_n = 0;
					w_selD_tmp = 1;
				end
				4'h7:	begin
					sel_mul_A_n = 10;
					sel_mul_B_n = 3;
					w_selD_tmp = 1;
				end
				4'h8:	begin
					sel_mul_A_n = 4;
					sel_mul_B_n = 4;
					w_selF_tmp = 1;
					sel_add_A_n = 3;
					sel_add_B_n = 3;
					w_selA_tmp = 2;
				end
				4'h9:	begin
					sel_mul_A_n = 1;
					sel_mul_B_n = 1;
					w_selB_tmp = 1;
					sel_add_A_n = 5;
					sel_add_B_n = 0;
					w_selA_tmp = 2;
					is_sub_n = 1;
				end
				4'hA:	begin
					sel_mul_A_n = 11;
					sel_mul_B_n = 1;
					w_selB_tmp = 1;
					sel_add_A_n = 3;
					sel_add_B_n = 0;
					w_selD_tmp = 2;
					is_sub_n = 1;
				end
				4'hB:	begin
					sel_mul_A_n = 3;
					sel_mul_B_n = 4;
					w_selD_tmp = 1;
				end
				4'hC:	begin
					sel_mul_A_n = 2;
					sel_mul_B_n = 2;
					w_selD_tmp = 1;
					sel_add_A_n = 3;
					sel_add_B_n = 1;
					w_selB_tmp = 2;
					is_sub_n = 1;
				end
			endcase
		PT_ADD:
			case(step_cnt)
				4'h0:	begin	// initialize 
					w_selG_tmp = 3;
					w_selH_tmp = 1;
					w_selI_tmp = 3;
				end
				4'h1:	begin
					sel_mul_A_n = 6;
					sel_mul_B_n = 3;
					w_selE_tmp = 1;
				end
				4'h2:	begin
					sel_mul_A_n = 3;
					sel_mul_B_n = 2;
					w_selD_tmp = 1;
				end
				4'h3:	begin
					sel_mul_A_n = 7;
					sel_mul_B_n = 3;
					w_selD_tmp = 1;
					sel_add_A_n = 0;
					sel_add_B_n = 4;
					w_selI_tmp = 2;
					is_sub_n = 1;
				end
				4'h4:	begin
					sel_add_A_n = 1;
					sel_add_B_n = 3;
					w_selF_tmp = 2;
					is_sub_n = 1;
				end
				4'h5:	begin
					sel_mul_A_n = 5;
					sel_mul_B_n = 5;
					w_selG_tmp = 1;
				end
				4'h6:	begin
					sel_mul_A_n = 8;
					sel_mul_B_n = 6;
					w_selE_tmp = 1;
					sel_add_A_n = 0;
					sel_add_B_n = 4;
					w_selA_tmp = 2;
				end
				4'h7:	begin
					sel_mul_A_n = 0;
					sel_mul_B_n = 4;
					w_selB_tmp = 1;
					sel_add_A_n = 3;
					sel_add_B_n = 1;
					w_selD_tmp = 2;
				end
				4'h8:	begin
					sel_mul_A_n = 8;
					sel_mul_B_n = 4;
					w_selE_tmp = 1;
					sel_add_A_n = 6;
					sel_add_B_n = 1;
					w_selA_tmp = 2;
					is_sub_n = 1;
				end
				4'h9:	begin
					sel_mul_A_n = 4;
					sel_mul_B_n = 3;
					w_selG_tmp = 1;
					sel_add_A_n = 0;
					sel_add_B_n = 0;
					w_selD_tmp = 2;
				end
				4'hA:	begin
					sel_mul_A_n = 2;
					sel_mul_B_n = 6;
					w_selC_tmp = 1;
					sel_add_A_n = 1;
					sel_add_B_n = 3;
					w_selB_tmp = 2;
					is_sub_n = 1;
				end
				4'hB:	begin
					sel_mul_A_n = 5;
					sel_mul_B_n = 1;
					w_selF_tmp = 1;
				end
				4'hC:	begin
					sel_mul_A_n = 0;
					sel_mul_B_n = 0;
					w_selD_tmp = 1;
					sel_add_A_n = 5;
					sel_add_B_n = 6;
					w_selF_tmp = 2;
					is_sub_n = 1;
				end
				4'hD:	begin
					sel_mul_A_n = 13;
					sel_mul_B_n = 5;
					w_selB_tmp = 1;
				end
			endcase
		CALXSQ:	begin
			sel_add_A_n = 0;
			sel_add_B_n = 0;
			w_selD_tmp = 1;
		end
		MULINV:
			case(step_cnt)
				4'd0:	begin
					w_selD_tmp = 3;		// D = MulyInv.r
					w_selC_tmp = 3;		// C = 1
				end
				1, 2, 3, 4, 5, 6, 7, 8, 9:	begin
					w_selC_tmp = 1;
					if (inv2power_mode == 0) begin
						sel_mul_A_n = 2;	// C = 2^(-k')
						sel_mul_B_n = 2;	
					end else begin
						sel_mul_A_n = 13;	// 2^(-1)
						sel_mul_B_n = 2;
					end
				end
				4'hA:	begin			// r * adjust
					sel_mul_A_n = 3;	// D = r
					sel_mul_B_n = 2;	// C = 2^(-k)
					w_selD_tmp = 1;		// D = Z^(-1)
				end
				4'hB:	begin			// p - mod(r * adjust)
					sel_add_A_n = 7;	// p
					sel_add_B_n = 3;	// D = mod(r * adjust)
					w_selD_tmp = 2;		// D = Z^(-1)
					is_sub_n = 1;
				end
				4'hC:	begin
					sel_mul_A_n = 3;	// D = Z^(-1)
					sel_mul_B_n = 3;
					w_selE_tmp = 1;		// E = Z^(-2)
				end
				4'hD:	begin
					sel_mul_A_n = 3;	// D = Z^(-1)
					sel_mul_B_n = 4;	// E = Z^(-2)
					w_selF_tmp = 1;		// F = Z^(-3)
				end
				4'hE:	begin
					sel_mul_A_n = 0;	// A = X0'
					sel_mul_B_n = 4;	// E = Z^(-2)
					w_selA_tmp = 1;		// A = X0
				end
				4'hF:	begin
					sel_mul_A_n = 1;	// B = Y0'
					sel_mul_B_n = 5;	// E = Z^(-3)
					w_selB_tmp = 1;		// B = Y0
				end
			endcase
		default: begin
			sel_add_A_n = 0;
			sel_add_B_n = 0;
			sel_mul_A_n = 0;
			sel_mul_B_n = 0;		
			is_sub_n = 0;
			w_selA_tmp = 0;
			w_selB_tmp = 0;
			w_selC_tmp = 0;
			w_selD_tmp = 0;
			w_selE_tmp = 0;
			w_selF_tmp = 0;
			w_selG_tmp = 0;
			w_selH_tmp = 0;
			w_selI_tmp = 0;
		end
	endcase
end

always @(posedge clk) begin	
	sel_add_A <= sel_add_A_n;
	sel_add_B <= sel_add_B_n;
	sel_mul_A <= sel_mul_A_n;
	sel_mul_B <= sel_mul_B_n;
	is_sub <= is_sub_n;
end

// Represent finished a calculation stage
always @(*) begin
	done_calc = 0;
	case (state)
		INIT:	done_calc = 1;
		PT_DBL: done_calc = valid_mod_prod;
		PT_ADD: 
			if (step_cnt == 0)
				done_calc = 1;
			else if (step_cnt == 4)
				done_calc = valid_mod_sum;
			else
				done_calc = valid_mod_prod;
		CALXSQ:
			done_calc = valid_mod_prod;
		MULINV:
			case(step_cnt)
				4'h0:	done_calc = valid_multinv;
				1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15:	
					done_calc = valid_mod_prod;	
				11:
					done_calc = valid_mod_sum;
			endcase
		default:
			done_calc = 0;
	endcase
end

// Register File Write Control ------------------------------------------------
always @(*)
	if (~done_calc) begin
		w_selA = 0;
		w_selB = 0;
		w_selC = 0;
		w_selD = 0;
		w_selE = 0;
		w_selF = 0;
		w_selG = 0;
		w_selH = 0;
		w_selI = 0;
	end else begin
		w_selA = w_selA_tmp;
		w_selB = w_selB_tmp;
		w_selC = w_selC_tmp;
		w_selD = w_selD_tmp;
		w_selE = w_selE_tmp;
		w_selF = w_selF_tmp;
		w_selG = w_selG_tmp;
		w_selH = w_selH_tmp;
		w_selI = w_selI_tmp;
	end
	
// Multiplicative Inverse -----------------------------------------------------
MultInv u_MultInv(
	.clk(clk),
	.en(en_multinv),
	.a(data[C]),
	.value(z_multinv),
	.power(powerk),
	.valid(valid_multinv)
);

endmodule