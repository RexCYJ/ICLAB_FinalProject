//======================================================================
//	Note:			Testbench for ECDH elliptic curve scalar
//	Author:			YuJia, Chen
//======================================================================

`timescale 1ns/100ps

`define TESTNUM 10

module test_batch;

parameter CYCLE = 10;
parameter MAXPATNUM = 100000;
parameter MAXCYCLE = 1000000000;
parameter BW_GF = 256;

integer i, j;

integer cycle_cnt;
integer start_cycle;
integer finish_cycle;
integer elaped_cycle;
integer system_start_cycle;
integer system_finish_cycle;

integer err_num;

// module ports
reg clk, srst_n;
reg start;
reg [BW_GF-1:0] k;
reg [BW_GF-1:0] Px;
reg [BW_GF-1:0] Py;
wire [BW_GF-1:0] Qx;
wire [BW_GF-1:0] Qy;
wire valid;

// input data file
reg [BW_GF-1:0] f_basepoint[0:2*MAXPATNUM-1];
reg [BW_GF-1:0] f_scalar[0:MAXPATNUM-1];
reg [BW_GF-1:0] f_endpoint[0:2*MAXPATNUM-1];

ECpoint_scalar u_top(
	.clk(clk),
	.rst_n(srst_n),
	.start(start),
	.k(k),
	.Px(Px),
	.Py(Py),
	.Qx(Qx),
	.Qy(Qy),
	.valid(valid)
);

reg is_show_detail;
reg is_show_simple;

initial begin
	`ifdef SHOWDETAIL
		is_show_detail = 1;
	`else
		is_show_detail = 0;
	`endif
	
	`ifdef SHOWSIMPLE
		is_show_simple = 1;
	`else
		is_show_simple = 0;
	`endif
end

always #(CYCLE/2) clk = ~clk;

// signal initialization
initial begin
	$display("\n\n\n>> TESTBENCH START ****************************");
	clk = 0;
	srst_n = 1;
	start = 0;
	$display(">> Reading: base point");
	$readmemh("./pat256/basepoint.dat", f_basepoint);
	$readmemh("./pat256/endpoint.dat", f_endpoint);
	$readmemh("./pat256/scalar.dat", f_scalar);
	
	@(negedge clk) srst_n = 0;
	@(negedge clk) srst_n = 1;
	@(negedge clk);

end

// feeding input signals
initial begin
	wait(~srst_n);
	wait(srst_n);
	
	$display("\n\n\n>> CRYPTO START ****************************");
	system_start_cycle = cycle_cnt;
	for (i = 0; i < `TESTNUM; i = i + 1) begin
		@(negedge clk);
		start = 1;
		start_cycle = cycle_cnt;
		Px = f_basepoint[0 + i*2];
		Py = f_basepoint[1 + i*2];
		k = f_scalar[i];
		@(negedge clk);
		start = 0;
		wait(valid);
		finish_cycle = cycle_cnt;
	end
	system_finish_cycle = cycle_cnt;
end

// cycle counter
initial begin
	cycle_cnt = 0;
	while(1) begin
		@(posedge clk);
		cycle_cnt = cycle_cnt + 1;
	end
end

reg err_flag;
real avg_cycle;

initial begin
	err_num = 0;
	for (j = 0; j < `TESTNUM; j = j + 1) begin
		err_flag = 0;
		wait(valid);
		#1;
		elaped_cycle = finish_cycle - start_cycle;

		if (Qx !== f_endpoint[0 + j*2] | Qy !== f_endpoint[1 + j*2]) begin
			err_num = err_num + 1;	
			err_flag = 1;
		end
		
		if (is_show_detail) begin
			$display("\n\n>> ================== Pattern No.%03d ==================", j);
			if (err_flag)
				$display(">> WRONG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
			$display(">> Golden: ");
			$display(">> \tX %64h", f_endpoint[0 + j*2]);
			$display(">> \tY %64h", f_endpoint[1 + j*2]);
			$display(">> RTL: ");
			$display(">> \tX %64h", Qx);
			$display(">> \tY %64h", Qy);
			$display(">> Total cycle:  %d", elaped_cycle);
		end
		if (is_show_simple) begin
			if (err_flag == 1) 
				$display(">> PAT No.%03d:  Wrong!", j);
			else
				$display(">> PAT No.%03d:  Correct!", j);
		end
		@(negedge clk);
		@(negedge clk);
	end
	
	avg_cycle = (system_finish_cycle - system_start_cycle) * 1.0 /(`TESTNUM);
	$display("\n\n>> Avg. Cycle:  %.2f", avg_cycle);

	if (err_num !== 0) begin
		$display("\n>> ============================================================");
		$display(">>      				%d	ERRORs!", err_num);
		$display(">> ============================================================\n");
	end else begin
		$display("\n>> ========!!VERIFICATION PASS, CONGRATULATION!!===============\n");
		$display("    ██████╗ ██████╗ ██████╗ ██████╗ ███████╗ ██████╗████████╗");
		$display("   ██╔════╝██╔═══██╗██╔══██╗██╔══██╗██╔════╝██╔════╝╚══██╔══╝");
		$display("   ██║     ██║   ██║██████╔╝██████╔╝█████╗  ██║        ██║   ");
		$display("   ██║     ██║   ██║██╔══██╗██╔══██╗██╔══╝  ██║        ██║   ");
		$display("   ╚██████╗╚██████╔╝██║  ██║██║  ██║███████╗╚██████╗   ██║   ");
		$display("    ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝ ╚═════╝   ╚═╝   ");
		$display(">> ===========================================================\n");
	end
	$finish;
end

// export waveform
// initial begin
// 	$fsdbDumpfile("batch.fsdb");
// 	$fsdbDumpvars("+mda");
// end

// terminate simulation if it takes too long
initial begin
	#(CYCLE * MAXCYCLE);
	$display("\n===================================================");
	$display("      Error!!! Simulation time is too long...      ");
	$display("   There might be something wrong in your code.    ");
 	$display("===================================================\n");
 	$finish;
end

endmodule