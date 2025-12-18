use std::error::Error;

/*
 Gaurav Sablok
 codeprog@icloud.com
*/

pub fn gcontent(inputstring: &str) -> Result<f32, Box<dyn Error>> {
    let inputstringchars = inputstring.chars().collect::<Vec<char>>();
    let mut finalgc: Vec<f32> = Vec::new();
    for i in inputstringchars.iter() {
        let mut count_a = 0usize;
        let mut count_t = 0usize;
        let mut count_g = 0usize;
        let mut count_c = 0usize;
        match i {
            'A' => count_a += 1usize,
            'T' => count_t += 1usize,
            'G' => count_g += 1usize,
            'C' => count_c += 1usize,
            _ => continue,
        }
        finalgc.push((count_g + count_c) as f32 / (count_a + count_t + count_g + count_c) as f32);
    }
    Ok(finalgc[0] as f32)
}
