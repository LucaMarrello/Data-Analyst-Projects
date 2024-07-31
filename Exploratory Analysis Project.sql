-- Exploratory Analysis

SELECT * 
FROM layoffs_staging2;

-- Exploring total laid off and percentage laid off variables to see maximum and minimum
-- number of people laid off in the dataset.

SELECT MAX(total_laid_off), MAX(percentage_laid_off)
FROM layoffs_staging2;

SELECT * 
FROM layoffs_staging2
WHERE percentage_laid_off = 1
ORDER BY funds_raised_millions DESC;

SELECT company, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company
ORDER BY 2 DESC;

SELECT MIN(`date`), MAX(`date`)
FROM layoffs_staging2;

SELECT industry, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY industry
ORDER BY 2 DESC;

SELECT country, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY country
ORDER BY 2 DESC;

SELECT YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY YEAR(`date`)
ORDER BY 1 DESC;

SELECT stage, SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY stage
ORDER BY 2 DESC;

-- Computing rolling total for people laid off by each month per year.

WITH Rolling_Total AS
(
SELECT substring(`date`, 1, 7) AS `Month`, SUM(total_laid_off) AS total_off
FROM layoffs_staging2
WHERE substring(`date`, 1, 7) IS NOT NULL
GROUP BY `MONTH`
ORDER BY 1 ASC
)
SELECT `Month`, total_off,
 SUM(total_off) OVER(ORDER BY `MONTH`) AS rolling_total
FROM Rolling_Total;

-- Exploring the companies who laid off the most people per year.

WITH company_Year (company, years, total_laid_off) AS
(
SELECT company, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY company, YEAR(`date`)
), Company_Year_Rank AS
(SELECT *, DENSE_RANK() OVER (PARTITION BY years ORDER BY total_laid_off DESC) AS
ranks
FROM company_Year
WHERE years IS NOT NULL
)
SELECT *
FROM Company_Year_Rank
WHERE ranks <= 5
;

-- Exploring the industries who laid off the most people by year.


WITH industry_Year (industry, years, total_laid_off) AS
(
SELECT industry, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY industry, YEAR(`date`)
), Industry_Year_Rank AS
(SELECT *, DENSE_RANK() OVER (PARTITION BY years ORDER BY total_laid_off DESC) AS
ranks
FROM industry_Year
WHERE years IS NOT NULL
)
SELECT *
FROM Industry_Year_Rank
WHERE ranks <= 5
;

-- Exploring the countries who laid off the most people by year.

WITH country_Year (country, years, total_laid_off) AS
(
SELECT country, YEAR(`date`), SUM(total_laid_off)
FROM layoffs_staging2
GROUP BY country, YEAR(`date`)
), Country_Year_Rank AS
(SELECT *, DENSE_RANK() OVER (PARTITION BY years ORDER BY total_laid_off DESC) AS
ranks
FROM country_Year
WHERE years IS NOT NULL
)
SELECT *
FROM Country_Year_Rank
WHERE ranks <= 5
;




