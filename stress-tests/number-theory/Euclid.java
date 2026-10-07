// Stress test for content/number-theory/Euclid.java. Not picked up by
// `make test` (which only runs *.cpp); run it from the repository root with
//   java stress-tests/number-theory/Euclid.java
// It compiles the header as it is in content/ and compares it to BigInteger.gcd.
import java.math.BigInteger;
import java.lang.reflect.Method;
import java.nio.file.*;
import java.util.*;

public class Euclid {
	static Method m;
	static BigInteger[] euclid(BigInteger a, BigInteger b) throws Exception {
		try {
			return (BigInteger[]) m.invoke(null, a, b);
		} catch (java.lang.reflect.InvocationTargetException e) {
			throw new RuntimeException("euclid(" + a + ", " + b + ") threw " + e.getCause());
		}
	}
	static void check(BigInteger a, BigInteger b) throws Exception {
		BigInteger[] r = euclid(a, b);
		String in = "euclid(" + a + ", " + b + ") = {" + r[0] + ", " + r[1] + ", " + r[2] + "}";
		if (!a.multiply(r[0]).add(b.multiply(r[1])).equals(r[2]))
			throw new RuntimeException("ax + by != d: " + in);
		if (!r[2].abs().equals(a.gcd(b)))
			throw new RuntimeException("|d| != gcd: " + in);
		if (a.signum() >= 0 && b.signum() >= 0 && r[2].signum() < 0)
			throw new RuntimeException("d < 0: " + in);
		// minimal coefficients
		if (a.signum() != 0 && b.signum() != 0 && !a.abs().equals(b.abs()))
			if (r[0].abs().multiply(r[2].abs()).compareTo(b.abs()) > 0 ||
					r[1].abs().multiply(r[2].abs()).compareTo(a.abs()) > 0)
				throw new RuntimeException("coefficients too large: " + in);
	}
	public static void main(String[] args) throws Exception {
		Path root = Paths.get(args.length > 0 ? args[0] : ".");
		String body = new String(Files.readAllBytes(
			root.resolve("content/number-theory/Euclid.java")));
		Path dir = Files.createTempDirectory("euclid");
		Files.write(dir.resolve("EuclidHeader.java"),
			("import java.math.BigInteger;\npublic class EuclidHeader {\n" + body + "}\n").getBytes());
		if (javax.tools.ToolProvider.getSystemJavaCompiler().run(null, null, null,
				"-d", dir.toString(), dir.resolve("EuclidHeader.java").toString()) != 0)
			throw new RuntimeException("header does not compile");
		Class<?> c = new java.net.URLClassLoader(new java.net.URL[]{dir.toUri().toURL()})
			.loadClass("EuclidHeader");
		m = c.getDeclaredMethod("euclid", BigInteger.class, BigInteger.class);
		m.setAccessible(true);

		Random rnd = new Random(12345);
		// exhaustive on small non-negative inputs
		for (int a = 0; a <= 300; a++) for (int b = 0; b <= 300; b++)
			check(BigInteger.valueOf(a), BigInteger.valueOf(b));
		// random non-negative inputs of every size
		for (int it = 0; it < 200000; it++) {
			BigInteger a = new BigInteger(rnd.nextInt(200), rnd), b = new BigInteger(rnd.nextInt(200), rnd);
			if (it % 3 == 0) { // force a large common factor
				BigInteger g = new BigInteger(rnd.nextInt(100), rnd);
				a = a.multiply(g); b = b.multiply(g);
			}
			check(a, b);
		}
		// worst case: consecutive Fibonacci numbers (~14400 steps), and big inputs
		BigInteger f0 = BigInteger.ZERO, f1 = BigInteger.ONE;
		for (int i = 0; i < 14400; i++) { BigInteger t = f0.add(f1); f0 = f1; f1 = t; }
		long st = System.nanoTime();
		check(f1, f0); check(f0, f1);
		for (int it = 0; it < 20; it++)
			check(new BigInteger(10000, rnd), new BigInteger(10000, rnd));
		System.out.println("large inputs: " + (System.nanoTime() - st) / 1000000 + " ms");
		// negative inputs: ax + by = d must still hold, with |d| = gcd(a, b)
		for (int a = -60; a <= 60; a++) for (int b = -60; b <= 60; b++)
			check(BigInteger.valueOf(a), BigInteger.valueOf(b));
		for (int it = 0; it < 100000; it++) {
			BigInteger a = new BigInteger(rnd.nextInt(200), rnd), b = new BigInteger(rnd.nextInt(200), rnd);
			if (rnd.nextBoolean()) a = a.negate();
			if (rnd.nextBoolean()) b = b.negate();
			check(a, b);
		}
		System.out.println("Tests passed!");
	}
}
