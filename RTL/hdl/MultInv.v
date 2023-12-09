//==================================================================================================
//  Note:          2022.12.31 Fall2022 NTHU ICLAB Final Project
//  Project:	   AES encryption and decryption with ECC key encryption
//  Author:	   	   Boris
//==================================================================================================

module MultInv #(
	parameter PAD = 2
)
(
	input clk,
	input srst_n,
	input en,
	input [`BW_GF-1:0] a,
	output reg [`BW_GF-1:0] value,
	output reg [9-1:0] power,
	output reg valid
);

integer i;

localparam BW = `BW_GF + PAD;

reg [`BW_GF-1:0] u, u_n;
reg [`BW_GF-1:0] v, v_n;
reg [`BW_GF+1-1:0] r, r_n;
reg [`BW_GF-1:0] s, s_n;
reg [9-1:0] k, k_n;
reg [`BW_GF-1:0] x, x_n;
// reg [`BW_GF-1:0] y, y_n;
reg busy, busy_n;
reg valid_n;
reg [`BW_GF-1:0] value_n;
reg [9-1:0] power_n;

reg signed [`BW_GF+1-1:0] x_pre, x_pre_n;
reg signed [`BW_GF+1-1:0] y_pre, y_pre_n;

wire signed [`BW_GF+1-1:0] p;

assign p = 257'h1_00000000_FFFFFFFE_FFFFFFFF_FFFFFFFF_FFFFFFFF_00000000_00000000_00000001;

always @(*) begin
	x_pre_n = {1'b1, u_n} + {1'b0, v_n};
	y_pre_n = r_n + s_n;
end

always @(posedge clk) begin
	x_pre <= x_pre_n;
	y_pre <= y_pre_n;
end

always @* begin
	u_n = u;
	v_n = v;
	r_n = r;
	s_n = s;
	k_n = k;
	x_n = x;
	// y_n = y;	// y doesn't be used
	power_n = power;
	value_n = value;
	
	valid_n = 0;

	if (busy) begin
		if (x != 0) begin
			busy_n = busy;
			if (u[0] == 0) begin
				u_n = {1'b1, u[`BW_GF-1:1]};	// u >>> 1;
				s_n = {s[`BW_GF-2:0], 1'b0};	// s <<< 1;
			end else if (v[0] == 0) begin
				v_n = {1'b0, v[`BW_GF-1:1]};	// v >>> 1;
				r_n = {r[`BW_GF-1:0], 1'b0}; 	// r <<< 1;
			end else begin
				x_n = x_pre[`BW_GF-1:0];
				// y_n = y_pre[`BW_GF-1:0];
				if (x_pre[`BW_GF] == 1) begin
					u_n = {x_pre[`BW_GF:1]};		// x_pre >>> 1
					r_n = {y_pre[`BW_GF:0]};		// y_pre
					s_n = {s[`BW_GF-2:0], 1'b0}; 	// s <<< 1;
				end else begin
					v_n = {x_pre[`BW_GF:1]};		// x_pre >>> 1;
					s_n = {y_pre[`BW_GF-1:0]};		// y_pre
					r_n = {r[`BW_GF-1:0], 1'b0};	// r <<< 1;
				end
			end
			k_n = k + 1;
		end else begin	// done
			busy_n = 0;
			valid_n = 1;
			power_n = k;
			if (r > `PRIME) value_n = r - `PRIME;
			else value_n = r;
		end
	end else begin
		busy_n = en;
		u_n = p[`BW_GF-1:0];
		v_n = a;
		r_n = 0;
		s_n = 1;
		k_n = 0;
		x_n = a;
		// y_n = 0;
		valid_n = 0;
		power_n = power;
		value_n = value;
	end
end

always @(posedge clk) begin
	u <= u_n;
	v <= v_n;
	r <= r_n;
	s <= s_n;
	k <= k_n;
	x <= x_n;
	// y <= y_n;
	// busy <= busy_n;
	// valid <= valid_n;
	value <= value_n;
	power <= power_n;
end

always @(posedge clk) begin
	if (~srst_n) begin
		valid <= 0;
		busy <= 0;
	end else begin
		valid <= valid_n;
		busy <= busy_n;
	end
end

endmodule
