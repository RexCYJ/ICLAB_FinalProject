//==================================================================================================
//  Note:			2022.12.31 Fall2022 NTHU ICLAB Final Project
//  Project:		AES encryption and decryption with ECC key encryption
//  Author:			Boris
//	Cycles:			2
//==================================================================================================

module mod_256 #(
	parameter a_len = 2*`BW_GF,
	parameter b_len = `BW_GF
)
(
	input clk,
	input [a_len-1:0] a,
	input valid,
	output reg [b_len-1:0] b,
	output reg finish
);

integer i;

localparam BITEXT = 4;
localparam LEN = `BW_GF+BITEXT;

reg [32-1:0] A [0:15];
wire signed [LEN-1-1:0] T;
wire signed [LEN-1-1:0] S1, S2, S3, S4;
wire signed [LEN-1-1:0] D1, D2, D3, D4;
reg signed [LEN-1-1:0] temp_b1, temp_b1_n;
reg signed [LEN-1-1:0] temp_b2, temp_b2_n;
reg signed [LEN-1:0] b1_adj, b2_adj;

reg signed [LEN-1:0] temp_b, temp_b_n;
reg signed [LEN-1:0] temp_b_q1, temp_b_q2;
reg signed [LEN-1:0] temp_b_q3, temp_b_q4;
reg signed [LEN-1:0] b_n;
reg valid_tmp1, valid_tmp2;

wire signed [LEN-1-1:0] p;

assign p = {3'b0, `PRIME};

always @(*)
	for (i = 0; i < 16; i = i + 1) 
		A[i] = a[i * 32 +: 32];

assign T  = {5'h0,  A[7],  A[6],  A[5],  A[4],  A[3],  A[2],  A[1],  A[0]};			// T

assign S1 = {4'h0, A[15], A[14], A[13], A[12], A[11], 32'h0, 32'h0, 32'h0, 1'b0};	// 2*S1
assign S2 = {4'h0, 32'h0, A[15], A[14], A[13], A[12], 32'h0, 32'h0, 32'h0, 1'b0};	// 2*S2
assign S3 = {5'h0, A[15], A[14], 32'h0, 32'h0, 32'h0, A[10], A[ 9], A[ 8]};			// S3
assign S4 = {5'h0, A[ 8], A[13], A[15], A[14], A[13], A[11], A[10], A[ 9]};			// S4

assign D1 = {5'h0, A[10], A[ 8], 32'h0, 32'h0, 32'h0, A[13], A[12], A[11]};	// 2*S2
assign D2 = {5'h0, A[11], A[ 9], 32'h0, 32'h0, A[15], A[14], A[13], A[12]};	// 2*S2
assign D3 = {5'h0, A[12], 32'h0, A[10], A[ 9], A[ 8], A[15], A[14], A[13]};	// 2*S2
assign D4 = {5'h0, A[13], 32'h0, A[11], A[10], A[ 9], 32'h0, A[15], A[14]};	// 2*S2

always @(*) begin
	temp_b1_n = T + S1 - D1 - D2 + p;
	temp_b2_n = S2 + S3 + S4 - D3 - D4 + p;
end

always @(posedge clk) begin
	temp_b1 <= temp_b1_n;
	temp_b2 <= temp_b2_n;
	valid_tmp1 <= valid;
end

always @(*) begin
	if (temp_b1 < 0)
		b1_adj = temp_b1 + p;
	else if (temp_b1 > p)
		b1_adj = temp_b1 - p;
	else
		b1_adj = temp_b1;
		
	if (temp_b2 < 0)
		b2_adj = temp_b2 + p;
	else if (temp_b2 > p)
		b2_adj = temp_b2 - p;
	else
		b2_adj = temp_b2;
end

always @(*)
	temp_b_n = b1_adj + b2_adj;

always @(posedge clk) begin
	temp_b <= temp_b_n;
	valid_tmp2 <= valid_tmp1;
end

always @(*) begin
	temp_b_q1 = temp_b - p;
	temp_b_q2 = temp_b - `PRIME * 2;	// -2p
	temp_b_q3 = temp_b - `PRIME * 3;	// -3p
	temp_b_q4 = temp_b - `PRIME * 4;	// -4p
end

always @(*) begin
	if (temp_b < {1'b0, p})
		b_n = temp_b;
	else if (temp_b_q1 < {1'b0, p})
		b_n = temp_b_q1;
	else if (temp_b_q2 < {1'b0, p})
		b_n = temp_b_q2;
	else if (temp_b_q3 < {1'b0, p})
		b_n = temp_b_q3;
	else
		b_n = temp_b_q4;
end

always @(posedge clk) begin
	b <= b_n;
	finish <= valid_tmp2;
end

endmodule
